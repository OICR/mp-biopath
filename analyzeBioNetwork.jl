#!/usr/bin/env julia

using JuMP
using Gurobi

m = Model(solver=GurobiSolver())
include("pi.pl")

UPPERBOUND = 200
NORMAL = UPPERBOUND/2
M = 1000

if length(ARGS) == 0
    println("Please provide the path to the pairwise interaction file in the first line")
    exit()
end

println("Performing analysis on pairwise interaction file: ", ARGS[1])

nodes = Pi.readFile(ARGS[1])

nodesList = collect(keys(nodes))

## Auxilary 
@variable(m, z[1:length(nodesList)], Bin)  #Above Normal
@variable(m, y[1:length(nodesList)], Bin)  #Below Normal
@variable(m, x[1:length(nodesList)] > 0, Int) ###actual nodes !!!

@variable(m, w[1:length(nodesList)] >= 0, Int)  #Above and optomized value
@variable(m, v[1:length(nodesList)] >= 0, Int)  #Below and optomized value

@variable(m, u[1:length(nodesList),1:UPPERBOUND], Bin)  #above being true 
@variable(m, t[1:length(nodesList),1:UPPERBOUND], Bin)  # neg being true

@variable(m, s[1:length(nodesList),1:UPPERBOUND,1:UPPERBOUND], Bin)

# Objective funtion
@objective(m, Min, sum{ M*(y[i] + z[i] + w[i] + v[i]) #auxilary variables
                       + y[i]*(NORMAL-x[i]) #weighting for above NORMAL
                       + z[i]*(x[i]-NORMAL) #weighting for below NORMAL
                        , i=1:length(nodesList)}
                       - sum{ u[i,j] + t[i,j], i=1:length(nodesList), 
                             i = 1:length(nodesList), j = 1:UPPERBOUND}
                       - sum{s[i,k,l],i=1:length(nodesList),
                             k = 1:UPPERBOUND, l = 1:UPPERBOUND})

@constraint(m, xzconst[i=1:length(nodesList)], z[i]*NORMAL >= x[i] - NORMAL)
@constraint(m, xyconst[i=1:length(nodesList)], y[i]*NORMAL >= NORMAL - x[i])

#=
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
=#


for nodeName in keys(nodes)
    currentIndexes = indexin([nodeName], nodesList)
    currentIndex = currentIndexes[1]

    ###Linearizing multiplication of and nodes
    relation = nodes[nodeName].relation
    parents = nodes[nodesName].parents
    if count(parents) == 0
        continue
    end

    parentIndexes = indexin([parents], nodesList)

    for idx = 1: UPPERBOUND
        @constraint(m,
                    ifabove[currentIndex,idx],
                    x[currentIndex] >= u[currentIndex,idx] * idx)
        @constraint(m,
                    ifbelow[currentIndex,idx],
                    idx >= x[currentIndex] * t[currentIndex,idx])
    end

    if count(parents) == 1
        @constraint(m,
                   andoneparentBelow[currentIndex],
                   x[parentIndexes[1]] + w[currentIndex] >= x[currentIndex])
        @constraint(m,
                   andoneparentAbove[currentIndex],
                   x[andPosParentIndexes[1]] - v[currentIndex] <= x[currentIndex])
    elseif length(orPosParents) > 0
        @constraint(m,
                    orposparent[currentIndex],
                    sum{x[orPosParentIndexes[a]], a=1:length(orPosParents)}/length(orPosParents)
                     + w[currentIndex] >= x[currentIndex])
        @constraint(m,
                    orposparent[currentIndex],
                    sum{x[orPosParentIndexes[a]], a=1:length(orPosParents)}/length(orPosParents)
                     - v[currentIndex] <= x[currentIndex])
    elseif numberofparents == length(andPosParents)
        for a = 1:UPPERBOUND, b = 1:UPPERBOUND 
            @constraint(m,
                        findalltrue[currentIndex,a,b],
                        u[andPosParentIndexes[1],a] + t[andPosParentIndexes[1],a]
                         + u[andPosParentIndexes[2],b] + t[andPosParentIndexes[2],b]
                         >= 4 * s[currentIndex,a,b])
        end
        @constraint(m,
                    setBasedOnParents[currentIndex],
                    sum{a*b*s[currentIndex,a,b], a=1:UPPERBOUND, b=1:UPPERBOUND}
                     + w[currentIndex] >= x[currentIndex])
        @constraint(m,
                    setBasedOnParents[currentIndex],
                    sum{a*b*s[currentIndex,a,b], a=1:UPPERBOUND, b=1:UPPERBOUND}
                     - v[currentIndex] <= x[currentIndex])
    end
end

println("solving model")
print(m)
solve(m)

for i in eachindex(nodesList)
        value = getvalue(x[i])
        println( i, "\t", nodesList[i], "\t\t", value)
end
exit()
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
