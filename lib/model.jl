using JuMP, AmplNLWriter
#model = Model(solver=CouenneNLSolver())
model = Model(solver=BonminNLSolver())

using JuMP
#using Gurobi
#model = Model(solver=GurobiSolver())


LB = 0.0001
UB = 10

weightRoot = 5
weightMeasured = 10000
weightHard = 10

if length(ARGS) == 0
    println("Please provide the path to the pairwise interaction file in the first line")
    exit()
end

println("Performing analysis on pairwise interaction file: ", ARGS[1])

nodes = Pi.readFile(ARGS[1])

nodesList = collect(keys(nodes))

measuredIdxs = Array{Integer}()

@variable(model, LB <= x[1:length(nodesList)] <= UB, start = 1)
@variable(model, LB <= x_bar[1:length(nodesList)] <= UB, start = 1)
@variable(model, LB <= p[1:length(nodesList)] <= UB)
@variable(model, LB <= n[1:length(nodesList)] <= UB)
@variable(model, m[1:length(measuredIdxs)])

rootIdxs = []
variableIdxs = []
for nodeName in keys(nodes)
    nodeIdxs = indexin([nodeName], nodesList)
    nodeIndex = nodeIdxs[1]

    measured = false
    j = 0
    for measuredIdx in measuredIdxs
        j = j + 1
        if measuredIdx == nodeIndex
            measured = true
            @constraint(model, measure[nodeIndex], m[j] == LB)
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
            posParentIdxs = indexin(nodes[nodeName].posParents, nodesList)
            @constraint(model,
                orposparentbelow[nodeIndex],
                sum{x[posParentIdxs[a]], a = 1:length(posParentIdxs)} / length(posParentIdxs) ==  x_bar[nodeIndex])
        end
    end
end

@objective(model,
           Min,
           weightHard * sum{p[variableIdxs[i]] + n[variableIdxs[i]], i = 1:length(variableIdxs)}
           + weightMeasured * sum{p[measuredIdxs[j]] + n[measuredIdxs[j]], j = 1:length(measuredIdxs)}
           + weightRoot * sum{p[rootIdxs[k]] + n[rootIdxs[k]], k = 1:length(rootIdxs)})

#println("solving model")
#print(model)
solve(model)

