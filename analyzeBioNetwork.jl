#!/usr/bin/env julia

include("pi.pl")

println("Running analysis on file:")

if length(ARGS) == 0
    println("Please provide the path to the pairwise interaction file in the first line")
else 
    println("Performing analysis on pairwise interaction file:")
    println(ARGS[1])

    pairwiseInteractions = Pi.readFile(ARGS[1])

    nodes = collection(keys(pairwiseInteractions))
    
 #   nodeList = Array{AbstractString,1}(node in nodes)

    println("node ", node)
println("type", typeof(node))
    i = 0
    for (nodeName, node) in pairwiseInteractions
        i += 0

        println("Node: ", getindex(nodeName, nodesList))
        if endswith(nodeName, "_RLE")
            println("\tIs a Reaction")
            if node.andOr == "AND"
                println("\t\t- AND relation")
                if length(node.negParents) > 0
                     println("\t\tHas a negative regulator ", length(node.negParents))
                     #Create constraint that flips negative regulator
                end
                #constraint 2 - take minimum from parents
            else 
                println("\t\t- OR relation")

                #contraint 3 
            end 
            #println(node.parents)
            #point

            #reaction

            ##positive
            ##negative

            #complex
            ##if all poisitive
        else
             println("is not a reaction")
             if node.andOr == "AND"
                println("\t\t- AND relation")
             else
                println("\t\t- OR relation")
             end
             #sets
              
             #positive
             #negative
             #contraint 6
             #contraint 7
             #contraint 8
             #contraint 9
             #contraint 11
             #contraint 12
        end
    end
end

