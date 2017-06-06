module NLmodel

using JuMP
using AmplNLWriter

function run(nodes, measurednodestate, keyoutputs, LB, UB, expression, verbose)
    model = Model(solver=CouenneNLSolver(["bonmin.nlp_log_level=0"; "bonmin.bb_log_level=0"]))

    weightRoot = 5
    weightMeasured = 10000
    weightHard = 10

    nodesList = collect(keys(nodes))

    @variable(model, LB <= x[1:length(nodesList)] <= UB, start = 1)
    @variable(model, LB <= x_bar[1:length(nodesList)] <= UB, start = 1)
    @variable(model, LB <= p[1:length(nodesList)] <= UB)
    @variable(model, LB <= n[1:length(nodesList)] <= UB)
    @variable(model, m[1:length(keys(measurednodestate))])

    rootIdxs = []
    variableIdxs = []
    for node in collect(keys(measurednodestate))
        if indexin([node],nodesList) == [0]
            delete!(measurednodestate, node)
        end
    end

    measuredIdxs = indexin(collect(keys(measurednodestate)), nodesList)

    for nodeName in keys(nodes)
        nodeIdxs = indexin([nodeName], nodesList)
        nodeIndex = nodeIdxs[1]

        measured = false
        j = 0
        for node in collect(keys(measurednodestate))
            j = j + 1
            if measuredIdxs[j] == nodeIndex
                measured = true
                rhs = measurednodestate[node] < LB? LB: measurednodestate[node]
                @constraint(model, measure[nodeIndex], m[j] == rhs)
                break
            end
        end

        if nodes[nodeName].relation == "ROOT"
            if measured
                @constraint(model,
                            pindexrootmeasured[nodeIndex],
                            p[nodeIndex] >= x[nodeIndex] - m[j])
                @constraint(model,
                            nindexrootmeasured[nodeIndex],
                            n[nodeIndex] >= m[j] - x[nodeIndex])
            else
                @constraint(model,
                            pindexroot[nodeIndex],
                            p[nodeIndex] >= x[nodeIndex] - 1)
                @constraint(model,
                            nindexroot[nodeIndex],
                            n[nodeIndex] >= 1 - x[nodeIndex])
                push!(rootIdxs, nodeIdxs[1])
            end
        else
            if measured
                @constraint(model,
                            pindexrootmeasured[nodeIndex],
                            p[nodeIndex] >= x[nodeIndex] - m[j])
                @constraint(model,
                            nindexrootmeasured[nodeIndex],
                            n[nodeIndex] >= m[j] - x[nodeIndex])
            else
                @constraint(model,
                            pindex[nodeIndex],
                            p[nodeIndex] >= x[nodeIndex] - x_bar[nodeIndex])
                @constraint(model,
                            nindex[nodeIndex],
                            n[nodeIndex] >= x_bar[nodeIndex] - x[nodeIndex])
            end

            push!(variableIdxs, nodeIdxs[1])

            if nodes[nodeName].relation == "AND"
                parentIndexes = indexin(nodes[nodeName].parents, nodesList)
                if length(parentIndexes) == 1
                    @constraint(model,
                                and[nodeIndex],
                                x[parentIndexes[1]] == x_bar[nodeIndex])
                else
                     @NLconstraint(model,
                                 and[nodeIndex],
                                 x[parentIndexes[1]] * x[parentIndexes[2]] == x_bar[nodeIndex])
                end
            elseif nodes[nodeName].relation == "ANDNEG"
                parentIndexes = indexin(nodes[nodeName].parents, nodesList)
                @NLconstraint(model,
                            and[nodeIndex],
                            x[parentIndexes[1]] / x[parentIndexes[2]] == x_bar[nodeIndex])
            elseif nodes[nodeName].relation == "OR"
                if verbose
                    println("Child node: $nodeName")
                end
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
                    @constraint(model,
                       orposparentbelow[nodeIndex],
                        sum{x[posParentIdxs[a]], a = 1:length(posParentIdxs)} / length(posParentIdxs) ==  x_bar[nodeIndex])
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

                    @constraint(model,
                        orposparentbelow[nodeIndex],
                        sum{evs[a] * x[posParentIdxs[a]], a = 1:length(posParentIdxs)} / total_expression ==  x_bar[nodeIndex])
                end
           end
        end
    end

    @NLobjective(model,
                 Min,
                 weightHard * sum{p[variableIdxs[i]] + n[variableIdxs[i]], i = 1:length(variableIdxs)}
                 + weightMeasured * sum{p[measuredIdxs[j]] + n[measuredIdxs[j]], j = 1:length(measuredIdxs)}
                 + weightRoot * sum{p[rootIdxs[k]] + n[rootIdxs[k]], k = 1:length(rootIdxs)})

    if verbose
        println("Solving Model")
        print(model)
    end

    status = solve(model)

    if verbose
        println("Objective value: ", getobjectivevalue(model))
    end

    keyresults = Dict()
    keyoutputsArray = collect(keyoutputs)
    if length(keyoutputs) > 0
        keyoutputIdx = indexin(keyoutputsArray, nodesList)
        for i in keyoutputIdx
            if i != 0
                keyresults[nodesList[i]] = getvalue(x[i])
            end
        end
    else
        for i in eachindex(nodesList)
            keyresults[nodesList[i]] = getvalue(x[i])
        end
    end

    x_values = Dict()
    x_bar_values = Dict()
    for i in eachindex(nodesList)
        x_values[nodesList[i]] = getvalue(x[i])
        x_bar_values[nodesList[i]] = getvalue(x_bar[i])
    end
 


    return [keyresults, x_values, x_bar_values]
end

end
