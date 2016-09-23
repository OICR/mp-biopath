#!/usr/bin/env julia

using JuMP
#using Gurobi
using NLopt
m = Model(solver=NLoptSolver(algorithm=:LD_SLSQP))
#m = Model(solver=GurobiSolver())
include("pi.pl")

println("Running analysis on file:")

UPPERBOUND = 6
NORMAL = UPPERBOUND/2
M = 1000

if length(ARGS) == 0
    println("Please provide the path to the pairwise interaction file in the first line")
    exit()
end

println("Performing analysis on pairwise interaction file: ", ARGS[1])

pairwiseInteractions = Pi.readFile(ARGS[1])
nodesList = collect(keys(pairwiseInteractions))



## Node
@variable(m, x[1:length(nodesList)] >=0, Int)

## Auxilary 
@variable(m, z[1:length(nodesList)], Bin)  #Above Normal
@variable(m, y[1:length(nodesList)], Bin)  #Below Normal
@variable(m, r[1:length(nodesList)] >= 0, Int)  #Or pos nodes

@variable(m, e[1:length(nodesList)] >= 0, Int)  #Above and value
@variable(m, f[1:length(nodesList)] >= 0, Int)  #Below and value

@variable(m, n[1:length(nodesList),1:UPPERBOUND], Bin)  #above being true 
@variable(m, l[1:length(nodesList),1:UPPERBOUND], Bin)  # neg being true

# Objective funtion
@objective(m, Min, sum{ M*(y[i] + z[i] + e[i] + f[i]) #auxilary variables
                       + y[i]*(NORMAL-x[i]) #weighting for above NORMAL
                       + z[i]*(x[i]-NORMAL) #weighting for below NORMAL
                        , i=1:length(nodesList)}
                       - sum{ n[i,j] + l[i,j], i=1:length(nodesList),
                             i=1:length(nodesList),j=1:UPPERBOUND})

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


for (nodeName, node) in pairwiseInteractions
    currentIndexes = indexin([nodeName], nodesList)
    currentIndex = currentIndexes[1]

    orPosParents =  collect(node.orPosParents)
    orPosParentIndexes = indexin(orPosParents, nodesList) 

    andPosParents =  collect(node.andPosParents)
    andPosParentIndexes = indexin(andPosParents, nodesList) 

    andNegParents =  collect(node.andNegParents)
    andNegParentIndexes = indexin(andNegParents, nodesList) 

    ###Linearizing multiplication of and nodes

    numberofparents = length(andPosParents) + length(andNegParents)

    if length(orPosParents) >= 1
        numberofparents = numberofparents + 1
    end

    if numberofparents == 0
        continue
    end

    idx = 0
    while idx < UPPERBOUND
        idx += 1
        @constraint(m,
                    ifabove[currentIndex,idx],
                    x[currentIndex] >= n[currentIndex,idx] * idx)
        @constraint(m,
                    ifbelow[currentIndex,idx],
                    idx >= x[currentIndex] * l[currentIndex,idx])
    end

    if numberofparents == 1
       @constraint(m,
                    andoneparentBelow[currentIndex],
                    x[andPosParentIndexes[1]] + e[currentIndex] >= x[currentIndex])
       @constraint(m,
                    andoneparentAbove[currentIndex],
                    x[andPosParentIndexes[1]] - f[currentIndex] <= x[currentIndex])
    elseif numberofparents == length(andPosParents)
       if numberofparents == 2 
            @NLconstraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND}
                + e[currentIndex] >= x[currentIndex])
            @NLconstraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND} 
                - f[currentIndex] <= x[currentIndex])
 #=       elseif numberofParents == 3
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND}
                + e[currentIndex] >= x[currentIndex])
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentsIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND} 
                - f[currentIndex] <= x[currentIndex])
       elseif numberofParents == 4
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND}
                + e[currentIndex] >= x[currentIndex])
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentsIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND} 
                - f[currentIndex] <= x[currentIndex])
       elseif numberofParents == 5
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND}
                + e[currentIndex] >= x[currentIndex])
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentsIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND} 
                - f[currentIndex] <= x[currentIndex])
       elseif numberofParents == 6
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[6],i]*l[andPosParentIndexes[6],i]/NORMAL*i, i=1:UPPERBOUND}
                + e[currentIndex] >= x[currentIndex])
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentsIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[6],i]*l[andPosParentIndexes[6],i]/NORMAL*i, i=1:UPPERBOUND} 
                - f[currentIndex] <= x[currentIndex])
       elseif numberofParents == 7
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[6],i]*l[andPosParentIndexes[6],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[7],i]*l[andPosParentIndexes[7],i]/NORMAL*i, i=1:UPPERBOUND}
                + e[currentIndex] >= x[currentIndex])
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentsIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[6],i]*l[andPosParentIndexes[6],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[7],i]*l[andPosParentIndexes[7],i]/NORMAL*i, i=1:UPPERBOUND} 
                - f[currentIndex] <= x[currentIndex])
       elseif numberofParents == 8
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[6],i]*l[andPosParentIndexes[6],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[7],i]*l[andPosParentIndexes[7],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[8],i]*l[andPosParentIndexes[8],i]/NORMAL*i, i=1:UPPERBOUND}
                + e[currentIndex] >= x[currentIndex])
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentsIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[6],i]*l[andPosParentIndexes[6],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[7],i]*l[andPosParentIndexes[7],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[8],i]*l[andPosParentIndexes[8],i]/NORMAL*i, i=1:UPPERBOUND} 
                - f[currentIndex] <= x[currentIndex])
       elseif numberofParents == 9
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[6],i]*l[andPosParentIndexes[6],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[7],i]*l[andPosParentIndexes[7],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[8],i]*l[andPosParentIndexes[8],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[9],i]*l[andPosParentIndexes[9],i]/NORMAL*i, i=1:UPPERBOUND}
                + e[currentIndex] >= x[currentIndex])
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentsIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[6],i]*l[andPosParentIndexes[6],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[7],i]*l[andPosParentIndexes[7],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[8],i]*l[andPosParentIndexes[8],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[9],i]*l[andPosParentIndexes[9],i]/NORMAL*i, i=1:UPPERBOUND} 
                - f[currentIndex] <= x[currentIndex])
       elseif numberofParents == 10
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[6],i]*l[andPosParentIndexes[6],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[7],i]*l[andPosParentIndexes[7],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[8],i]*l[andPosParentIndexes[8],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[9],i]*l[andPosParentIndexes[9],i]/NORMAL*i, i=1:UPPERBOUND}
                *sum{n[andPosParentIndexes[10],i]*l[andPosParentIndexes[10],i]/NORMAL*i, i=1:UPPERBOUND}
                + e[currentIndex] >= x[currentIndex])
            @constraint(m,
                setBasedOnParents[currentIndex],
                sum{n[andPosParentIndexes[1],i]*l[andPosParentsIndexes[1],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[2],i]*l[andPosParentIndexes[2],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[3],i]*l[andPosParentIndexes[3],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[4],i]*l[andPosParentIndexes[4],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[5],i]*l[andPosParentIndexes[5],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[6],i]*l[andPosParentIndexes[6],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[7],i]*l[andPosParentIndexes[7],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[8],i]*l[andPosParentIndexes[8],i]/NORMAL*i, i=1:UPPERBOUND} 
                *sum{n[andPosParentIndexes[10],i]*l[andPosParentIndexes[10],i]/NORMAL*i, i=1:UPPERBOUND} 
                - f[currentIndex] <= x[currentIndex])
=#
        end
    elseif length(orPosParents) > 0
  #=       @constraint(m, 
                    orposreg[currentIndex], 
                    sum{x[parentIndex],
                    parentIndex=orPosParentIndexes}/length(orPosParentIndexes)  == r[currentIndex])
        if length(orPosParents) == numberofparents
            @constraint(m,
                        onlyOrPosParentBelow[currentIndex],
                        r[currentIndex] + e[currentIndex] >= x[currentIndex])
            @constraint(m,
                        onlyOrPosparentAbove[currentIndex],
                        r[currentIndex] - f[currentIndex] <= x[currentIndex])
        elseif length(andPosParents) > 0 && length(andNegParents) == 0
        else #andPos and orNeg an orPos
        end
    elseif length(orPosParents) >= 0 && length(andNegParents) >= 0
=#
    else #this is for nodes with no parents
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
