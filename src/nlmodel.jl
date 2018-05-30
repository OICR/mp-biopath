module NLmodel

using JuMP
using PiecewiseLinearOpt
#using Gurobi
using AmplNLWriter
using CoinOptServices

function runModel(nodes, measuredNodeStateFull, UB, expression, options, verbose)
    model = Model(solver=AmplNLSolver(CoinOptServices.couenne, options))
#    model = Model(solver=GurobiSolver(Presolve=0))
#    model = Model(solver = OsilSolver(solver = "couenne", options))
    cnt = 0
    LB = 0.01
    weightRoot = 5
    weightMeasured = 10000
    weightHard = 10

    nodesList = collect(keys(nodes))
    measuredNodeState = filter((node,v)->indexin([node], nodesList) != [0], measuredNodeStateFull)
    zDiscreteNodes = vcat(-UB:1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1:UB)
    yDiscreteNodes = vcat(0.1 , 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1:(2*UB))
    discreteValues = [0.01,2,100] #[0, 0.39, 0.78, 1.56, 3.125, 6.25, 12.5, 25, 50, UB] #vcat(1:0.1:1,2:7:UB)

    @variables model begin
        x[1:length(nodesList)], (lowerbound = LB, start = 1, upperbound = UB)
        x_bar[1:length(nodesList)], (lowerbound = LB, start = 1, upperbound = UB)
        p[1:length(nodesList)], (lowerbound = 0, upperbound = UB)
        n[1:length(nodesList)], (lowerbound = 0, upperbound = UB)
        z[1:length(nodesList)], (lowerbound = 0, upperbound = UB)
 #       λ[1:length(nodesList),1:length(zDiscreteNodes)], (lowerbound = 0, upperbound = 1)
 #       w[1:length(nodesList),1:length(zDiscreteNodes)], Bin
 #       y[1:length(zDiscreteNodes)], Bin
 #       Z[1:length(nodesList),1:length(zDiscreteNodes)], (lowerbound = 0, start = 0, upperbound = 2 * UB)
 #       Y[1:length(nodesList),1:length(yDiscreteNodes)], (lowerbound = -UB, start = 0, upperbound = UB)
 #       YR[1:length(nodesList),1:length(yDiscreteNodes)], Bin
 #       ZR[1:length(nodesList),1:length(zDiscreteNodes)], Bin
 #       ySquared[1:length(nodesList),1:length(yDiscreteNodes)], (lowerbound = 0, start = 0, upperbound = 2 * UB)
 #       zSquared[1:length(nodesList),1:length(zDiscreteNodes)], (lowerbound = 0, start = 0, upperbound = 2 * UB)
    end

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
                rhs = measuredNodeState[node] < LB? LB: measuredNodeState[node]
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
                push!(rootIdxs, nodeIdxs[1])
            end
        else
            if measured
                @constraint(model, p[nodeIndex] >= x[nodeIndex] - m[j])
                @constraint(model, n[nodeIndex] >= m[j] - x[nodeIndex])
            else
                @constraint(model, p[nodeIndex] >= x[nodeIndex] - x_bar[nodeIndex])
                @constraint(model, n[nodeIndex] >= x_bar[nodeIndex] - x[nodeIndex])
            end

            push!(variableIdxs, nodeIdxs[1])

            if nodes[nodeName].relation == "AND"
                parentIndexes = indexin(nodes[nodeName].parents, nodesList)
                if length(parentIndexes) == 1
                    @constraint(model, x[parentIndexes[1]] == x_bar[nodeIndex])
                else
                    # Piecewise linearize x[parentIndexes[1]] * x[parentIndexes[2]] = x_bar[nodeIndex]

      #              @constraint(model, x[parentIndexes[1]] - x[parentIndexes[2]] == Z[nodeIndex])
      #              @constraint(model, x[parentIndexes[1]] + x[parentIndexes[2]] == Y[nodeIndex])
 
                    f(u,v) = u * v
                    println("printing")
                    println(x[parentIndexes[1]])
                    println(x[parentIndexes[2]])
                    println(discreteValues)
                    println(f)
                    z[nodeIndex] = piecewiselinear(model, x[parentIndexes[1]], x[parentIndexes[2]], discreteValues, discreteValues, f, method=:ZigZag)
                    @constraint(model, z[nodeIndex] == x_bar[nodeIndex])

#                    @NLconstraint(model, x[parentIndexes[1]] * x[parentIndexes[2]] == x_bar[nodeIndex])
#                    lastz = 0.00
#                    println("typeof")
#                    println(typeof(lastz))
#	            for zIndex in eachindex(zDiscreteNodes)
#                        z = zDiscreteNodes[zIndex]
#                        slope = (z ^ 2 - lastz ^ 2) / z - lastz
#			b = z ^ 2 - (slope * z)
#                        @constraint(model, slope * Z[nodeIndex] + b == zSquared[nodeIndex, zIndex])
#                        @constraint(model, ZR[nodeIndex, zIndex] * (lastz + 0.00001) <= Z[nodeIndex] <= ZR[nodeIndex, zIndex] * z)
#                        lastz = z
#        	    end
#            
#		    lasty = - UB - 1
#                    for yIndex in eachindex(yDiscreteNodes)
#                        y = yDiscreteNodes[yIndex]
#
#                        slope = (y ^ 2 - lasty ^ 2) / y - lasty
#			b = y ^ 2 - (slope * y)
#
#                        @constraint(model, slope * Y[nodeIndex, yIndex] + b == ySquared[nodeIndex, yIndex])
#                        @constraint(model, YR[nodeIndex,yIndex] * (lasty + 0.00001) <= Y[nodeIndex, yIndex] <= YR[nodeIndex, yIndex] * yIndex)
#
#                        lasty = y
#		    end
# 	
#                   @constraint(model, sum(YR[nodeIndex,i] for i in 2:length(yDiscreteNodes)) == 1)
#	           @constraint(model, sum(ZR[nodeIndex,i] for i in 1:length(zDiscreteNodes)) == 1)
#        
#                   @constraint(model, (sum(YR[nodeIndex,i] * ySquared[nodeIndex,i] for i in 1:length(yDiscreteNodes)) 
#					 - sum(ZR[nodeIndex,i] * zSquared[nodeIndex,i] for i in 1:length(zDiscreteNodes))) / 4 == x_bar[nodeIndex])
                end
            elseif nodes[nodeName].relation == "NEG"
                parentIndexes = indexin(nodes[nodeName].parents, nodesList)
                @NLconstraint(model, 1 / x[parentIndexes[1]] == x_bar[nodeIndex])
                #    fdivideone(u) = 1 / u
                #    z = piecewiselinear(model, x[parentIndexes[1]], discreteValues, fdivideone)
                #    @constraint(model, z == x_bar[nodeIndex])
            elseif nodes[nodeName].relation == "ANDNEG"
                parentIndexes = indexin(nodes[nodeName].parents, nodesList)
               #  fdivide(u,v) = u / v
               #  z = piecewiselinear(model, x[parentIndexes[1]], x[parentIndexes[2]], discreteValues, discreteValues, fdivide)
               #  @constraint(model, z == x_bar[nodeIndex])
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

    @objective(model,
                 Min,
                 weightHard * sum(p[variableIdxs[i]] + n[variableIdxs[i]] for i = 1:length(variableIdxs))
                 + weightMeasured * sum(p[measuredIdxs[j]] + n[measuredIdxs[j]] for j = 1:length(measuredIdxs))
                 + weightRoot * sum(p[rootIdxs[k]] + n[rootIdxs[k]] for k = 1:length(rootIdxs)))

    if verbose
        println("Solving Model")
        print(model)
    end

    status = solve(model)

    println("status")
    println(status)
    if verbose
        println("Objective value: ", getobjectivevalue(model))
    end

    results = Dict()
    for i in eachindex(nodesList)
         results[nodesList[i]] = getvalue(x[i])
    end

    x_values = Dict()
    x_bar_values = Dict()
    for i in eachindex(nodesList)
        x_values[nodesList[i]] = getvalue(x[i])
        x_bar_values[nodesList[i]] = getvalue(x_bar[i])
    end
    zs = getvalue(z[1])
    zss = getvalue(z[2])
    zsss = getvalue(z[3])
    println(z)
    println("z: $zs $zss $zsss")
    println("results: $results")
    println("x_values: $x_values")
    println("x_bar_values: $x_bar_values")
    return [results, x_values, x_bar_values]
end

end
