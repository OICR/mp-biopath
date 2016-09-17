#!/usr/bin/env julia

using JuMP
using Gurobi

include("pi.pl")

println("Running analysis on file:")

NORMAL = 100
M = 1000

if length(ARGS) == 0
    println("Please provide the path to the pairwise interaction file in the first line")
    exit()
end

println("Performing analysis on pairwise interaction file:")
println(ARGS[1])

pairwiseInteractions = Pi.readFile(ARGS[1])

nodes = keys(pairwiseInteractions)
   
#nodeList is created in order to find and index for a node by name. This is used in order to make 
#variable names e.g. X1 (X[i])

nodesList = Array{AbstractString, 1}()
for node in nodes
    push!(nodesList, node)

    for parent in pairwiseInteractions[node].negParents
        indexin([parent], nodesList) == 0 || push!(nodesList, parent)
    end

    for parent in pairwiseInteractions[node].posParents
        indexin([parent], nodesList) == 0 || push!(nodesList, parent)
    end
end


# Preparing an optimization model
m = Model(solver=GurobiSolver())

# Declaring variables

##node variables
@variable(m, x[1:length(nodesList)] >=0, Int)

#these are auxilary for deterining if x is above or below NORMAL
@variable(m, xz[1:length(nodesList)], Bin)
@variable(m, xy[1:length(nodesList)], Bin)


###Objective funtion
@objective(m, Min, sum{ M*(xy[i]+xz[i]) #auxilary variables
                       + xy[i]*(NORMAL-x[i]) #weighting for above NORMAL
                       + xz[i]*(x[i]-NORMAL) #weighting for below NORMAL
                        , i=1:length(nodesList)})

@constraint(m, xzconst[i=1:length(nodesList)], xz[i]*NORMAL >= x[i] - NORMAL)
@constraint(m, xyconst[i=1:length(nodesList)], xy[i]*NORMAL >= NORMAL - x[i])


#Signalling by ERBB2
@constraint(m, fixed, x[71] == 0)  #downregulating EGFR and ERBB3
#expecting 73 and 73 to be down

for (nodeName, node) in pairwiseInteractions
    currentIndexes = indexin([nodeName], nodesList)
    currentIndex = currentIndexes[1]

    if node.andOr == "AND"
        #println("AND relation")
        if length(node.negParents) > 0
            #println("-- Has a negative regulator(s): ", length(node.negParents))
            #Create constraint that flips negative regulator
        else
            println(node.posParents)
            for posParent in node.posParents
                parentIndexes = indexin([posParent], nodesList)
                println(typeof(parentIndexes))
                parentIndex = parentIndexes[1]
                @constraint(m, andposregcurrentIndex[parentIndex], x[currentIndex] <= x[parentIndex])

            end
        end
    else 
        #println("OR relation")
        if length(node.negParents) > 0
             #println("-- Has a negative regulator(s): ", length(node.negParents))
             #Create constraint that flips negative regulator
        else
             #println("-- All positive regulator(s): ", length(node.posParents))
        end

        #contraint 3 
    end 
end

solve(m)
print(m)
println("Optimal Solutions:")
for i in eachindex(nodesList)
    println( i, ") ", nodesList[i], " = ", getvalue(x[i]))
end
