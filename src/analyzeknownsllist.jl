module AnalyzeKnownSLList

include("../src/pi.jl")
include("../src/evidence.jl")
include("../src/nlmodel.jl")
include("../src/keyoutputs.jl")
include("../src/essential.jl")
include("../src/results.jl")
include("../src/sl.jl")
include("../src/idMap.jl")
#include("../src/coexpress.jl")


function run(pifile, slfilename, lowerbound, upperbound, dbidfile, verbose)
    pinodes = Pi.readFile(pifile)
    slnodes = SL.getNodes(dbidfile)
    pinodesSet = Set(keys(pinodes))

#    coexpress = CoExpress.getForPiNodes(pinodesSet, dbidfile)
#    println(coexpress)
#    exit()

    essentialgenes = Essential.getGenes(collect(keys(pinodes)), dbidfile)
    nodetogene     = DbIdNameMapping.nodeToGene(dbidfile)

    slessential = Dict{String,Any}[]
    slpinodes = String[]
    slpinodesPairs = []
    for sl in slnodes
        for nodea in sl["GeneANodes"]
            for nodeb in sl["GeneBNodes"]
                if in(nodea, pinodesSet) && in(nodeb, pinodesSet) && pinodes[nodea].relation == "ROOT" && pinodes[nodeb].relation == "ROOT"
                    if in(nodeb, slpinodes) == false
                        push!(slpinodes, nodeb)
                    end
                    if in(nodea, slpinodes) == false
                        push!(slpinodes, nodea)
                    end
                    push!(slpinodesPairs, [nodea, nodeb])
                    for essentialgene in essentialgenes
                        if nodetogene[essentialgene] == nodetogene[nodea] || nodetogene[essentialgene] == nodetogene[nodeb]
                            continue
                        end

                        slcopy = copy(sl)
                        delete!(slcopy, "GeneANodes")
                        delete!(slcopy, "GeneBNodes")

                        slcopy["GeneANode"] = nodea
                        slcopy["GeneBNode"] = nodeb
                        slcopy["EssentialGene"] = essentialgene
                        push!(slessential, slcopy)
                    end
                end
            end
        end
    end
    allGenes = DbIdNameMapping.allGeneReferenceProduct(dbidfile)
    allGenesSet = Set(allGenes)
    for node in slpinodes
        (sampleresults, x, x_bar) = NLmodel.run(pinodes,
                                                Dict(node => 0),
                                                essentialgenes,
                                                lowerbound,
                                                upperbound,
                                                verbose)

        for (essentialnode, essentialvalue) in sampleresults
            for sl in slessential
                if sl["EssentialGene"] == essentialnode
                    if sl["GeneANode"] == node
                        sl["GeneANodeValue"] = essentialvalue
                    elseif sl["GeneBNode"] == node
                        sl["GeneBNodeValue"] = essentialvalue
                    end
                end
            end
        end
    end

    sloutfile = open(slfilename, "w")
    headercolumns = ["Count"]
    for column in keys(slessential[1])
        push!(headercolumns, column)
    end
    header = join(headercolumns, "\t")
    write(sloutfile, string(header, "\n"))
    flush(sloutfile)

    for slnodespair in slpinodesPairs
        (sampleresults, x, x_bar) = NLmodel.run(pinodes,
                                                Dict(slnodespair[1] => 0,
                                                     slnodespair[2] => 0),
                                                     essentialgenes,
                                                     lowerbound,
                                                     upperbound,
                                                     verbose)

         for (essentialnode, essentialvalue) in sampleresults
             for sl in slessential
                 if sl["EssentialGene"] == essentialnode && sl["GeneANode"] == slnodespair[1] && sl["GeneBNode"] == slnodespair[2]
                     sl["EssentialGeneValue"] = essentialvalue
                 end
             end
         end
    end

    sloutfile = open(slfilename, "w")
    headercolumns = ["Count"]
    for column in keys(slessential[1])
        push!(headercolumns, column)
    end
    header = join(headercolumns, "\t")
    write(sloutfile, string(header, "\n"))

    count = 0
    for sl in slessential
        count = count + 1
        columnvalues = String[]
        for column in headercolumns
            if column == "Count"
                push!(columnvalues, string(count))
            else
                push!(columnvalues, string(sl[column]))
            end
        end
        write(sloutfile, join(columnvalues, "\t"), "\n")
    end

    close(sloutfile)
end

end
