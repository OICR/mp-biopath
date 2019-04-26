module NLmodel

using JuMP
#using AmplNLWriter, 
using Ipopt
#using CoinOptServices

function runModel(nodes, measuredNodeStateFull, LB, UB, expression, options, verbose)
    model = Model(with_optimizer(Ipopt.Optimizer))

    weightRoot = 5
    weightMeasured = 10000
    weightHard = 100

    nodesList = collect(keys(nodes))
    measuredNodeState = filter((args)->indexin([args[1]], nodesList) != [0], measuredNodeStateFull)

    @variable(model, LB <= x[1:length(nodesList)] <= UB, start = 1)
    @variable(model, LB <= x_bar[1:length(nodesList)] <= UB, start = 1)
    @variable(model, LB <= p[1:length(nodesList)] <= UB)
    @variable(model, LB <= n[1:length(nodesList)] <= UB)
    
    numberMeasuredNodes = length(keys(measuredNodeState))

    if numberMeasuredNodes > 0
       @variable(model, m[1:numberMeasuredNodes])
    end
 
    rootIdxs = []
    variableIdxs = []
    measuredIdxs = indexin(collect(keys(measuredNodeState)), nodesList)
    for nodeName in keys(nodes)
        nodeIdxs = indexin([nodeName], nodesList)
        nodeIndex = nodeIdxs[1]

        measured = false
        j = 0
        for node in collect(keys(measuredNodeState))
            j = j + 1
            if measuredIdxs[j] == nodeIndex && nodes[nodeName].relation == "ROOT"
                measured = true
                rhs = measuredNodeState[node] < LB ? LB : measuredNodeState[node]
                @constraint(model, m[j] == rhs)
                break
            end
        end

        if nodes[nodeName].relation == "ROOT"
            if measured
                @constraint(model, p[nodeIndex] >= x[nodeIndex] - m[j])
                @constraint(model, n[nodeIndex] >= m[j] - x[nodeIndex])
            else
                @constraint(model, p[nodeIndex] >= x[nodeIndex] - 1)
                @constraint(model, n[nodeIndex] >= 1 - x[nodeIndex])
                push!(rootIdxs, nodeIndex)
            end
        else
            if measured
                @constraint(model, p[nodeIndex] >= x[nodeIndex] - m[j])
                @constraint(model, n[nodeIndex] >= m[j] - x[nodeIndex])
            else
                @constraint(model, p[nodeIndex] >= x[nodeIndex] - x_bar[nodeIndex])
                @constraint(model, n[nodeIndex] >= x_bar[nodeIndex] - x[nodeIndex])
            end
            push!(variableIdxs, nodeIndex)

            if nodes[nodeName].relation == "AND"
                parentIndexes = indexin(nodes[nodeName].parents, nodesList)
                if length(parentIndexes) == 1
                    @constraint(model, x[parentIndexes[1]] == x_bar[nodeIndex])
                else
                    @NLconstraint(model, x[parentIndexes[1]] * x[parentIndexes[2]] == x_bar[nodeIndex])
                end
            elseif nodes[nodeName].relation == "NEG"
                parentIndexes = indexin(nodes[nodeName].parents, nodesList)
                @NLconstraint(model, 1 / x[parentIndexes[1]] == x_bar[nodeIndex])
            elseif nodes[nodeName].relation == "ANDNEG"
                parentIndexes = indexin(nodes[nodeName].parents, nodesList)
                @NLconstraint(model, x[parentIndexes[1]] / x[parentIndexes[2]] == x_bar[nodeIndex])
            elseif nodes[nodeName].relation == "OR"
                posParentIdxs = indexin(nodes[nodeName].posParents, nodesList)
                count_expression = 0
                total_expression = 0
                for parent in nodes[nodeName].posParents
                    if haskey(expression, parent) && expression[parent] != ""
                        total_expression += parse(Float64, expression[parent])
                        count_expression += 1
                        ev = expression[parent]
                        if verbose
                            println("Parent node: $parent\tExpression value: $ev")
                        end
                    end
                end
                if verbose
                    println("")
                end

                if count_expression == 0
                    @constraint(model, sum(x[posParentIdxs[a]] for a = 1:length(posParentIdxs)) / length(posParentIdxs) ==  x_bar[nodeIndex])
                else
                    average_expression = total_expression / count_expression
                    evs = []
                    for parent in nodes[nodeName].posParents
                        if haskey(expression, parent) && expression[parent] != ""
                            push!(evs, parse(Float64, expression[parent]))
                        else
                            push!(evs, average_expression);
                        end
                    end

                    total_expression = sum(evs)

                    @constraint(model, sum(evs[a] * x[posParentIdxs[a]] for a = 1:length(posParentIdxs)) / total_expression ==  x_bar[nodeIndex])
                end
           end
        end
    end

    if !isempty(measuredIdxs)
        measuredIdxs = measuredIdxs[measuredIdxs .!= nothing]
    end

    @NLobjective(model,
                 Min,
                 weightHard * sum(p[variableIdxs[i]] + n[variableIdxs[i]] for i = 1:length(variableIdxs))^2
                 + weightMeasured * sum(p[measuredIdxs[j]] + n[measuredIdxs[j]] for j = 1:length(measuredIdxs))^2
                 + weightRoot * sum(p[rootIdxs[k]] + n[rootIdxs[k]] for k = 1:length(rootIdxs))^2)

    if verbose
        println("Solving Model")
        println(model)
    end

    optimize!(model)

    if verbose
        println("Objective value: ", getobjectivevalue(model))
    end

    results = Dict()
    for i in eachindex(nodesList)
         results[nodesList[i]] = JuMP.value(x[i])
    end

    x_values = Dict()
    x_bar_values = Dict()
    for i in eachindex(nodesList)
        x_values[nodesList[i]] = JuMP.value(x[i])
        x_bar_values[nodesList[i]] = JuMP.value(x_bar[i])
    end

    return [results, x_values, x_bar_values]
end

end
