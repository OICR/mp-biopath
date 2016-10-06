#!/usr/bin/env julia

using JuMP

using JuMP, AmplNLWriter
#model = Model(solver=CouenneNLSolver())
model = Model(solver=BonminNLSolver())

#using Gurobi
#model = Model(solver=GurobiSolver())
include("../lib/pi.jl")

LB = 0.0001
UB = 10

weightRoot = 5
weightMeasured = 10000
weightHard = 1000

if length(ARGS) == 0
    println("Please provide the path to the pairwise interaction file in the first line")
    exit()
end

println("Performing analysis on pairwise interaction file: ", ARGS[1])

nodes = Pi.readFile(ARGS[1])

nodesList = collect(keys(nodes))


measuredIdxs = Array{Integer}()

if ismatch(r"test", ARGS[1])
    measuredIdxs = indexin(["b", "c", "a"], nodesList)
elseif ismatch(r"Signaling_by_ERBB2", ARGS[1])
    measuredIdxs = indexin(["1963571", "p-10Y-ERBB3-1_[plasma_membrane]_54424"], nodesList)
elseif ismatch(r"DNA_Double-Strand_Break_Repair", ARGS[1])
    measuredIdxs = indexin(["RAD52_[nucleoplasm]_62640",
#                    "PALB2_[nucleoplasm]_241569",
                   "BRCA2_[nucleoplasm]_50952"
                   ], nodesList)
elseif ismatch(r"PIP3_activates_AKT_signaling", ARGS[1])
    #test AKT upregulated
    #measuredidxs = indexin(["AKT1_[cytosol]_58253", "AKT2_[cytosol]_49860", "AKT3_[cytosol]_415917"], nodesList)

    #test PIK3CA upgregulated test
    #measuredidxs = indexin(["PIK3CA_[cytosol]_61074"], nodesList)

    #test3 PTEN downregulated
    measuredIdxs = indexin(["PTEN_mRNA_[cytosol]_2318745"], nodesList)
elseif ismatch(r"DNA_Damage_Reversal", ARGS[1])
    measuredIdxs = indexin(["5657646"], nodesList)
end

@variable(model, LB <= x[1:length(nodesList)] <= UB, start = 1)
@variable(model, LB <= x_bar[1:length(nodesList)] <= UB, start = 1)
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
        if measuredIdx == nodeIdxs
            measured = true
            @constraint(model, measure[nodeIndex], m[j] == 1)
            break
        end
    end

    if nodes[nodeName].relation == "ROOT"
        if !measured
           push!(rootIdxs, nodeIdxs[1])
        end
    else
        push!(variableIdxs, nodeIdxs[1])
        parentIndexes = indexin(nodes[nodeName].parents, nodesList)
        if nodes[nodeName].relation == "AND"
            if length(parentIndexes) == 1
                @constraint(model,
                            and[nodeIndex],
                            x[parentIndexes[1]] == x_bar[nodeIndex])
            else
                 @constraint(model,
                             and[nodeIndex],
                             x[parentIndexes[1]] * x[parentIndexes[2]] == x_bar[nodeIndex])
            end
        elseif nodes[nodeName].relation == "ANDNEG"
            @constraint(model,
                        and[nodeIndex],
                        x[parentIndexes[1]] / x[parentIndexes[2]] == x_bar[nodeIndex])
        elseif nodes[nodeName].relation == "OR"
            posParentIdxs = indexin(nodes[nodeName].posParents, nodesList)
            @constraint(m,
                orposparentbelow[currentIndex],
                sum{x[posParentIdxs[a]], a = 1:length(posParentIdxs)} / length(posParentIdxs) ==  x_bar[nodeIndex])
        end
    end
end

@objective(model,
           Min,
           sum{weightHard * (x[variableIdxs[i]] - x_bar[variableIdxs[i]])^2, i = 1:length(variableIdxs)}
           + sum{weightMeasured * (x[measuredIdxs[j]] - m[j])^2, j = 1:length(measuredIdxs)}
           + sum{weightRoot * (x[rootIdxs[k]] - 1)^2, k = 1:length(rootIdxs)})

println("solving model")
#print(model)
#exit()
solve(model)

function valueToState(value)
    if 1 > value
        return "Down Regulated"
    elseif 1 < value
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
