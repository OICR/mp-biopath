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

println("Performing analysis on pairwise interaction file: ", ARGS[1])

pairwiseInteractions = Pi.readFile(ARGS[1])
nodesList = collect(keys(pairwiseInteractions))

m = Model(solver=GurobiSolver())

# Variables

## Node
@variable(m, x[1:length(nodesList)] >=0, Int)

## Auxilary 
@variable(m, z[1:length(nodesList)], Bin)  #Above Normal
@variable(m, y[1:length(nodesList)], Bin)  #Below Normal
@variable(m, o[1:length(nodesList)], Int)  #Min or max value of the node without neg reg
@variable(m, r[1:length(nodesList)], Int)  #for or pos nodes
@variable(m, n[1:length(nodesList)], Int)  #for neg regulation

# Objective funtion
@objective(m, Min, sum{ M*(y[i]+z[i]) #auxilary variables
                       + y[i]*(NORMAL-x[i]) #weighting for above NORMAL
                       + z[i]*(x[i]-NORMAL) #weighting for below NORMAL
                        , i=1:length(nodesList)})

@constraint(m, xzconst[i=1:length(nodesList)], z[i]*NORMAL >= x[i] - NORMAL)
@constraint(m, xyconst[i=1:length(nodesList)], y[i]*NORMAL >= NORMAL - x[i])

if ismatch(r"Signaling_by_ERBB2", ARGS[1])
    idxs = indexin(["1963571", "p-10Y-ERBB3-1_[plasma_membrane]_54424"], nodesList)
    @constraint(m, EGFR, x[idxs[1]] == 0)
    @constraint(m, ERBB3, x[idxs[2]] == 0)
elseif ismatch(r"DNA_Double-Strand_Break_Repair", ARGS[1])
    idxs = indexin(["RAD52_[nucleoplasm]_62640",
#                    "PALB2_[nucleoplasm]_241569",
                   "BRCA2_[nucleoplasm]_50952"
                   ], nodesList)

    @constraint(m, RAD52, x[idxs[1]] == 0) 
    @constraint(m, secondone, x[idxs[2]] == 0) 
end

for (nodeName, node) in pairwiseInteractions
    currentIndexes = indexin([nodeName], nodesList)
    currentIndex = currentIndexes[1]

    orPosParents =  collect(node.orPosParents)
    orPosParentIndexes = indexin(orPosParents, nodesList) 
    andPosParents =  collect(node.andPosParents)
    andPosParentIndexes = indexin(andPosParents, nodesList) 

hello = "    for (i, idx) in enumerate(andPosParentIndexes)
        if i == 1
           product = product * x[idx]/100
        else
           product = (product) * (x[idx]/100)
        end
    end "

    if length(orPosParentIndexes) > 0
        @constraint(m, 
          orposreg[currentIndex], 
          sum{x[parentIndex],
              parentIndex=orPosParentIndexes}/length(orPosParentIndexes)  == r[currentIndex])
    else
         if length(andPosParentIndexes) >= 2
         @NLconstraint(m, andposreg[currentIndex], x[andPosParentIndexes[1]] >= x[currentIndex])
         end
    end

    andNegParents =  collect(node.andNegParents)
    andNegParentIndexes = indexin(andNegParents, nodesList) 
    if length(andNegParentIndexes) > 0
        @constraint(m, 
          andnegreg[currentIndex], 
          sum{x[parentIndex],
              parentIndex=andNegParentIndexes} - length(andNegParentIndexes)*NORMAL == n[currentIndex])
    else
        @constraint(m, n[currentIndex] == 100)
    end

end

solve(m)
#print(m)

function valueToState(value)
    if 70 > value
        return "Down Regulated"
    elseif 130 < value
        return "Up Regulated"
    else
        return "Normal"
    end
end

println("Optimal Solutions:")

if ismatch(r"Signaling_by_ERBB2", ARGS[1])
    println("EGFR path")
    for i in eachindex(nodesList)
        if indexin([nodesList[i]], ["1963571",
                                    "1963589_RLE",
                                    "1963573",
                                    "1963582_RLE",
                                    "1963585",
                                    "1963588",
                                    "1963578_RLE",
                                    "1248746",
                                    "1250195_RLE",
                                    "1250194",
                                    "1250486_RLE",
                                    "1250479",
                                    "1250463_RLE",
                                    "109783" #lethal node
                                   ]) != [0]
            value = getvalue(x[i])
            println( i, "\t", nodesList[i], "\t\t", value, "\t", valueToState(value))
        end
    end
    
    println()
    println("ERBB3 path")
    for i in eachindex(nodesList)
        if indexin([nodesList[i]], ["p-10Y-ERBB3-1_[plasma_membrane]_54424",
                                    "1248743",
                                    "1248749",
                                    "1963572",
                                    "1250189_RLE",
                                    "1250508",
                                    "1250462",
                                    "PI(3_4_5)P3_[plasma_membrane]_188411" #end result
                                    ]) != [0]
            value = getvalue(x[i])
            println( i, "\t", nodesList[i], "\t\t", value, "\t", valueToState(value))
        end
    end

elseif ismatch(r"DNA_Double-Strand_Break_Repair", ARGS[1]) 
   println("rad52")
   for i in eachindex(nodesList)
        if indexin([nodesList[i]], ["rad52_[nucleoplasm]_62640",
                                    "75998_RLE",
                                    "83899",
                                    "5686587_RLE",
                                    "5693580_RLE",
                                    "5693564_RLE",
                                    "5693590",
                                    "5686642_RLE",
                                    "5686613",
                                    "5686657_RLE",
                                    "5686662",
                                    "5686663_RLE",
                                    "5686663" #lethal node
                                    ]) != [0]
            value = getvalue(x[i])
            println( i, "\t", nodesList[i], "\t\t", value, "\t", valueToState(value))
        end
   end

   println()

   println("PALB2")
   for i in eachindex(nodesList)
        if indexin([nodesList[i]], ["PALB2_[nucleoplasm]_241569",
                                    "5693620_RLE",
                                    "5685838_RLE",
                                    "5685826",
                                    "5693593_RLE",
                                    "5686104",
                                    "5686440_RLE", #slit over these three reactions
                                    "5693539_RLE", 
                                    "5693589_RLE",
                                    "5686432",
                                    "5686469_RLE", #lethal node
                                    "5686228",
                                    "5693584_RLE",
                                    "5686493",
                                    "5686483_RLE", #lethal node
                                    "5686410_RLE", #lethal node
                                    "84009",
                                    "5693558_RLE", #lethal node
                                    ]) != [0]
            value = getvalue(x[i])
            println( i, "\t", nodesList[i], "\t\t", value, "\t", valueToState(value))
        end
   end

   println()

   println("BRCA2")
   for i in eachindex(nodesList)
        if indexin([nodesList[i]], ["BRCA2_[nucleoplasm]_50952",
                                    "5685242_RLE",
                                    "5693561_RLE",
                                    "d"
                                    ]) != [0]
            value = getvalue(x[i])
            println( i, "\t", nodesList[i], "\t\t", value, "\t", valueToState(value))
        end
    end
end
