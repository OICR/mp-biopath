#!/usr/bin/env julia

using JuMP
using Gurobi

include("pi.pl")

println("Running analysis on file:")

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
i = 0
nodesList = Array{AbstractString, 1}()
for node in nodes
    push!(nodesList, node)

    for parent in pairwiseInteractions[node].negParents
        indexin([parent], nodesList) != 0 || push!(nodesList, parent)
    end

    for parent in pairwiseInteractions[node].posParents
        indexin([parent], nodesList) != 0 || push!(nodesList, parent)
    end
end
   
for (nodeName, node) in pairwiseInteractions
    currentNodeIndex = indexin([nodeName], nodesList)

    if node.andOr == "AND"
        println("AND relation")
        if length(node.negParents) > 0
             #println("-- Has a negative regulator(s): ", length(node.negParents))
             #Create constraint that flips negative regulator
        else
             #println("-- All positive regulator(s): ", length(node.posParents))
             #constraint 2 - take minimum from parents
        end
    else 
        println("OR relation")
        if length(node.negParents) > 0
             #println("-- Has a negative regulator(s): ", length(node.negParents))
             #Create constraint that flips negative regulator
        else
             #println("-- All positive regulator(s): ", length(node.posParents))
        end

        #contraint 3 
    end 
end

