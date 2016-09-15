#!/usr/bin/env julia

include("pi.pl")

println("Running analysis on file:")

if length(ARGS) == 0
    println("Please provide the path to the pairwise interaction file in the first line")
else 
    println("Performing analysis on pairwise interaction file:")
    println(ARGS[1])

    pairwiseInteractions  = Pi.readFile(ARGS[1])

    nodes = keys(pairwiseInteractions)

    i = 0
    for (nodeName, node) in pairwiseInteractions
        i += 0
        if node.andOr == "AND"
            if endswith(nodeName, "_RLE")
                println("reaction: ", nodeName)
                #constraint 1 - if one down other up
                #constraint 2 - take minimum from parents
            else 
                println("Normal reaction: ", nodeName)

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
             if endswith(nodeName, "_RLE")
                 println("reaction: ", nodeName)
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

