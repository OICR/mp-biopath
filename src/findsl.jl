module FindSL

include("valuetostate.jl")
include("dbidnamemapping.jl")
include("essential.jl")
include("expression.jl")
include("nlmodel.jl")
include("pi.jl")

function run(pifile, slfilename, lowerbound, upperbound, downregulatedcutoff, upregulatedcutoff, pairwisefile, dbidfile, tissueType, verbose)
    pinodes = Pi.readFile(pifile)
    expression = Expression.get(dbidfile, tissueType)
    #essentialgenes = Essential.getGenes(collect(keys(pinodes)), dbidfile)
    essentialgenes = AbstractString["419195","5668934"]
    allGenes = DbIdNameMapping.allGeneReferenceProduct(dbidfile)
    allGenesSet = Set(allGenes)

    geneList = ["CIT_[cytosol]_53072",
        "DOCK1_[cytosol]_89242",
	"DOCK2_[cytosol]_89244",
	"ERBB2_[plasma_membrane]_54422",
	"HGF(32-494)_[extracellular_region]_56498",
	"HGF(495-728)_[extracellular_region]_56498",
	"KALRN_[cytosol]_401664",
	"MET_[plasma_membrane]_402617",
	"p-Y1234_Y1235_Y1349_Y1356-MET_[plasma_membrane]_402617",
	"p-Y1420-DCC_[plasma_membrane]_53436",
	"p-Y829-TIAM1_[cytosol]_195051",
	"TIAM1_[cytosol]_195051",
	"TRIO_[cytosol]_401916"]

    downregulatednodestates = []
    for node in keys(pinodes)
        if in(node, allGenesSet) == false
            continue
        end

        if pinodes[node].relation != "ROOT"
            continue
        end
        if contains(node, "PSEUDONODE") || in(node, Set(essentialgenes))
            continue
        end

        for i in [0,2]
            if in(node, Set(essentialgenes)) == false || i == 2
                (sampleresults, x, x_bar) = NLmodel.run(pinodes,
                                                        Dict(node => i),
                                                        essentialgenes,
                                                        lowerbound,
                                                        upperbound,
                                                        expression,
                                                        verbose)

                essentialnodes = []
                for resultnode in keys(sampleresults)
                    value = sampleresults[resultnode]
                    state = ValueToState.getState(value,
                                                  downregulatedcutoff,
                                                  upregulatedcutoff)

                    if state == "Down Regulated"
                        push!(essentialnodes, Dict("value" => value,
                                                   "name" => resultnode ))
                    end
                end
                if length(essentialnodes) != 0
                    nodestate = Dict("name" => node,
                                     "state" => i,
                                     "essentialNodes" => essentialnodes)
                    push!(downregulatednodestates, nodestate)
                end
            end
        end
    end

    sloutfile = open(slfilename, "w")
    write(sloutfile, "count\tnode_one\tstate\tnode_one_value\tnode_two\tstate\tnode_two_value\teffected_node\teffected_node_value\n")
    flush(sloutfile)

    count = 1
    candidates = 0
    index = 0
    indextwo = 0
    for nodeone in downregulatednodestates
        index = index + 1
        for nodetwo in downregulatednodestates
            indextwo = indextwo + 1
            if indextwo <= index
                continue
            end

            if in(nodeone["name"], geneList) == false && in(nodetwo["name"], geneList) == false
               continue
            end

            candidates = candidates + 1

            (sampleresults, x, x_bar) = NLmodel.run(pinodes,
                                            Dict(nodeone["name"] => nodeone["state"],
                                                 nodetwo["name"] => nodetwo["state"]),
                                                 essentialgenes,
                                                 lowerbound,
                                                 upperbound,
                                                 expression,
                                                verbose)
            for resultnode in keys(sampleresults)
               value = sampleresults[resultnode]
               state = ValueToState.getState(value,
                                        downregulatedcutoff,
                                        upregulatedcutoff)

               if state == "Down Regulated"
                   essentialnodesone = nodeone["essentialNodes"]
                   for essentialnodeone in essentialnodesone
                        if essentialnodeone["name"] != resultnode
                            continue
                        end
                        essentialnodestwo = nodetwo["essentialNodes"]
                        for essentialnodetwo in essentialnodestwo
                            if essentialnodetwo["name"] != resultnode
                                continue
                            end

                            nameone = essentialnodeone["name"]
                            valueone = essentialnodeone["value"]
                            nametwo = essentialnodetwo["name"]
                            valuetwo = essentialnodetwo["value"]

                            if (value < valueone) && (value < valuetwo)
                                write(sloutfile, "$count\t")
                                write(sloutfile, nodeone["name"])
                                write(sloutfile, "\t")
                                write(sloutfile, nodeone["state"])
                                write(sloutfile, "\t")
                                write(sloutfile, string(valueone))
                                write(sloutfile, "\t")
                                write(sloutfile, nodetwo["name"])
                                write(sloutfile, "\t")
                                write(sloutfile, nodetwo["state"])
                                write(sloutfile, "\t")
                                write(sloutfile, string(valuetwo))
                                write(sloutfile, "\t")
                                write(sloutfile, "$resultnode\t$value\n")
                                flush(sloutfile)

                            end
                            count = count + 1
                            break
                        end

                    end
                end
            end
        end
    end
    println("number of candidate SL $candidates\n")
    close(sloutfile)
end

end
