#!/usr/bin/env julia

using JuMP
using Gurobi

m = Model(solver=GurobiSolver())
include("../lib/pi.jl")

UPPERBOUND = 10
NORMAL = UPPERBOUND/2
M = 100000

if length(ARGS) == 0
    println("Please provide the path to the pairwise interaction file in the first line")
    exit()
end

println("Performing analysis on pairwise interaction file: ", ARGS[1])

nodes = Pi.readFile(ARGS[1])

nodesList = collect(keys(nodes))

## Auxilary 
@variable(m, z[1:length(nodesList)], Bin)                           #Above Normal
@variable(m, y[1:length(nodesList)], Bin)                           #Below Normal
@variable(m, 0 <= x[1:length(nodesList)] <= UPPERBOUND, Int, start = NORMAL)        ###actual nodes !!!
@variable(m, 0 <= w[1:length(nodesList)] <= UPPERBOUND, Int, start = 0)        #Above and optomized value
@variable(m, 0 <= v[1:length(nodesList)] <= UPPERBOUND, Int, start = 0)        #Below and optomized value
@variable(m, u[1:length(nodesList),1:UPPERBOUND], Bin)              #Above being true 
@variable(m, t[1:length(nodesList),1:UPPERBOUND], Bin)              #Neg being true
@variable(m, s[1:length(nodesList),1:UPPERBOUND,1:UPPERBOUND], Bin) #Binary of two parent values

@objective(m, Min, sum{ 8 * NORMAL * (y[i] + z[i]) + M * (w[i] + v[i]) #auxilary variables
                       + 8 * y[i]*(NORMAL - x[i]) #weighting for above NORMAL
                       + 8 * z[i]*(x[i] - NORMAL) #weighting for below NORMAL
                        , i=1:length(nodesList)}
                       - sum{ 8 * u[i,j], i = 1:length(nodesList), j = 1:UPPERBOUND}
                       - sum{ 8 * t[i,j], i = 1:length(nodesList), j = 1:UPPERBOUND}
                       - sum{s[i,k,l],i = 1:length(nodesList),
                             k = 1:UPPERBOUND, l = 1:UPPERBOUND})
# Observed
observedidxs = Array{Integer}()

if ismatch(r"test", ARGS[1])
    observedidxs = indexin(["b", "c", "a"], nodesList)
    @constraint(m, test, x[observedidxs[1]] == 10)
    @constraint(m, test, x[observedidxs[2]] == 5)
elseif ismatch(r"Signaling_by_ERBB2", ARGS[1])
    observedidxs = indexin(["1963571", "p-10Y-ERBB3-1_[plasma_membrane]_54424"], nodesList)
    @constraint(m, EGFR, x[observedidxs[1]] == 0)
    @constraint(m, ERBB3, x[observedidxs[2]] == 0)
elseif ismatch(r"DNA_Double-Strand_Break_Repair", ARGS[1])
    idxs = indexin(["RAD52_[nucleoplasm]_62640",
#                    "PALB2_[nucleoplasm]_241569",
                   "BRCA2_[nucleoplasm]_50952"
                   ], nodesList)

    @constraint(m, RAD52, x[observedidxs[1]] == 0) 
    @constraint(m, secondone, x[observedidxs[2]] == 0) 
elseif ismatch(r"PIP3_activates_AKT_signaling", ARGS[1])
 #test AKT upregulated
 #  observedidxs = indexin(["AKT1_[cytosol]_58253", "AKT2_[cytosol]_49860", "AKT3_[cytosol]_415917"], nodesList)
 #  @constraint(m, AKT1, x[observedidxs[1]] == UPPERBOUND)
 #  @constraint(m, AKT2, x[observedidxs[2]] == UPPERBOUND)
 #  @constraint(m, AKT3, x[observedidxs[2]] == UPPERBOUND)

 #test PIK3CA upgregulated test
 #  observedidxs = indexin(["PIK3CA_[cytosol]_61074"], nodesList)
 #  @constraint(m, PIK3CA, x[observedidxs[1]] == UPPERBOUND)

 #test3 PTEN downregulated
   observedidxs = indexin(["PTEN_mRNA_[cytosol]_2318745", "PTEN_Gene_[nucleoplasm]_5632940"], nodesList)
   @constraint(m, PTENcytosol, x[observedidxs[1]] == 1)
   @constraint(m, PTENnucleoplasm, x[observedidxs[2]] == 1)


elseif ismatch(r"DNA_Damage_Reversal", ARGS[1])
   println("Setting obser for DNA_Damage_Reversal")
   observedidxs = indexin(["5657646"], nodesList)
   @constraint(m, oneetadsDNA, x[observedidxs[1]] == 1)
end

for nodeName in keys(nodes)
    currentIndexes = indexin([nodeName], nodesList)
    currentIndex = currentIndexes[1]

    observedIndex = indexin([currentIndex], observedidxs)

    for idx = 1: UPPERBOUND
        @constraint(m,
            ifabove[currentIndex,idx],
            x[currentIndex] >= u[currentIndex,idx] * idx)
        @constraint(m,
            ifbelow[currentIndex,idx],
            idx >= x[currentIndex] * t[currentIndex,idx])
    end

    if nodes[nodeName].relation == "ROOT" && observedIndex[1] == 0
        @constraint(m, xzconst[currentIndex], z[currentIndex] * NORMAL >= x[currentIndex] - NORMAL)
        @constraint(m, xyconst[currentIndex], y[currentIndex] * NORMAL >= NORMAL - x[currentIndex])
        continue
    else
        @constraint(m, yzxconstzero[currentIndex], z[currentIndex] + y[currentIndex] == 0)
    end

    if (nodes[nodeName].relation != "OR" && (length( nodes[nodeName].parents) == 1))
        parentIndexes = indexin(nodes[nodeName].parents, nodesList)
        if observedIndex[1] != 0
            @constraint(m,
                andoneparentBelow[currentIndex],
                x[parentIndexes[1]] + w[currentIndex] >= x[currentIndex])
            @constraint(m,
                andoneparentAbove[currentIndex],
                x[parentIndexes[1]] - v[currentIndex] <= x[currentIndex])
        else
            @constraint(m, wvconstzero[currentIndex], x[currentIndex] == x[parentIndexes[1]])
            @constraint(m, wvconstzero[currentIndex], w[currentIndex] + v[currentIndex] == 0)
        end
    elseif nodes[nodeName].relation == "AND" || nodes[nodeName].relation == "ANDNEG"
        parentIndexes = indexin(nodes[nodeName].parents, nodesList)
        for a = 1:UPPERBOUND, b = 1:UPPERBOUND 
            @constraint(m,
                findalltrue[currentIndex,a,b],
                u[parentIndexes[1],a] + t[parentIndexes[1],a]
                 + u[parentIndexes[2],b] + t[parentIndexes[2],b]
                 >= 4 * s[currentIndex,a,b])
        end
           
        if true || observedIndex[1] != 0
            if nodes[nodeName].relation == "AND" 
                @constraint(m,
                    setBasedOnParents[currentIndex],
                    sum{(a / NORMAL) * (b / NORMAL) * s[currentIndex,a,b], a = 1:UPPERBOUND, b = 1:UPPERBOUND} * NORMAL
                             + w[currentIndex] >= x[currentIndex])
                @constraint(m,
                    setBasedOnParents[currentIndex],
                    sum{(a / NORMAL)*(b / NORMAL) * s[currentIndex,a,b], a = 1:UPPERBOUND, b = 1:UPPERBOUND} * NORMAL
                     - v[currentIndex] <= x[currentIndex])
            else #ANDNEG
                @constraint(m,
                    setBasedOnParents[currentIndex],
                    sum{(a / NORMAL) / (b / NORMAL) * s[currentIndex,a,b], a = 1:UPPERBOUND, b = 1:UPPERBOUND} * NORMAL
                     + w[currentIndex] >= x[currentIndex])
                @constraint(m,
                    setBasedOnParents[currentIndex],
                    sum{(a / NORMAL) / (b / NORMAL) * s[currentIndex,a,b], a = 1:UPPERBOUND, b = 1:UPPERBOUND} * NORMAL
                     - v[currentIndex] <= x[currentIndex])
            end
        else
             if nodes[nodeName].relation == "AND" 
                @constraint(m,
                    setBasedOnParents[currentIndex],
                    sum{(a / NORMAL) * (b / NORMAL) * s[currentIndex,a,b], a = 1:UPPERBOUND, b = 1:UPPERBOUND} * NORMAL  == x[currentIndex])
            else #ANDNEG
                @constraint(m,
                    setBasedOnParents[currentIndex],
                    sum{(a / NORMAL) / (b / NORMAL) * s[currentIndex,a,b], a = 1:UPPERBOUND, b = 1:UPPERBOUND} * NORMAL == x[currentIndex])
            end
        end
        @constraint(m,
                    totals[currentIndex],
                    sum{s[currentIndex,a,b], a = 1:UPPERBOUND, b = 1:UPPERBOUND} == 1)
    elseif nodes[nodeName].relation == "OR"
        posParentIdxs = indexin(nodes[nodeName].posParents, nodesList)
        negParentIdxs = indexin(nodes[nodeName].negParents, nodesList)
        if true || observedIndex[1] != 0
            @constraint(m,
                orposparentbelow[currentIndex],
                ((sum{x[posParentIdxs[a]], a = 1:length(posParentIdxs)} + sum{UPPERBOUND - x[negParentIdxs[b]], b = 1:length(negParentIdxs)}) / (length(posParentIdxs) + length(negParentIdxs)))
                 + w[currentIndex] >= x[currentIndex])
            @constraint(m,
                orposparentabove[currentIndex],
                ((sum{x[posParentIdxs[a]], a = 1:length(posParentIdxs)} + sum{UPPERBOUND - x[negParentIdxs[b]], b = 1:length(negParentIdxs)}) / (length(posParentIdxs) + length(negParentIdxs)))
                 - v[currentIndex] <= x[currentIndex])
        else
            @constraint(m,
                orposparentbelow[currentIndex],
                sum{x[parentIdxs[a]], a = 1:length(parents)} / length(parents) == x[currentIndex])
        end
    end
end

println("solving model")
#print(m)
#exit()
solve(m)

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
elseif ismatch(r"PIP5_activates_AKT", ARGS[1])
    println("PIP3 results")
    for i in eachindex(nodesList)
        value = getvalue(x[i])
        if value != NORMAL
            println( i, "\t", nodesList[i], "\t\t", value)
        end
    end
else
    println("else")
    for i in eachindex(nodesList)
        value = getvalue(x[i])
        println( i, "\t", nodesList[i], "\t\t", value, "\t", valueToState(value))
    end
#=    println(nodesList)
    println("x")
    println(getvalue(x))
    println("y")
    println(getvalue(y))
    println("z")
    println(getvalue(y))
    println("v")
    println(getvalue(v))
    println("w")
    println(getvalue(w))
    println("u")
    println(getvalue(u))
    println("t")
    println(getvalue(t))
=#

 #   println("s")
 #   println(getvalue(s))


end
