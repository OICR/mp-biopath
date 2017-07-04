module FindSL

include("valuetostate.jl")
include("dbidnamemapping.jl")
include("essential.jl")
include("nlmodel.jl")
include("pi.jl")

function run(pifile, slfilename, lowerbound, upperbound, downregulatedcutoff, upregulatedcutoff, pairwisefile, dbidfile, tissueType, verbose)
    pinodes = Pi.readFile(pifile)
    expression = Expression.get(dbidfile, tissueType)
    essentialgenes = Essential.getGenes(collect(keys(pinodes)), dbidfile)
    allGenes = DbIdNameMapping.allGeneReferenceProduct(dbidfile)
    allGenesSet = Set(allGenes)

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
    write(sloutfile, "count\tnode_one\tnode_one_value\tnode_two\tnode_two_value\teffected_node\teffected_node_value\n")
    flush(sloutfile)

    count = 1
    candidates = 0
    index = 0
    indextwo = 0
    for nodeone in downregulatednodestates
        index = index + 1
        for nodetwo in downregulatednodestates
            found = false
            indextwo = indextwo + 1
            if indextwo <= index
                continue
            end
            candidates = candidates + 1

            sampleresults = NLmodel.run(pinodes,
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
                                write(sloutfile, string(valueone))
                                write(sloutfile, "\t")
                                write(sloutfile, nodetwo["name"])
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
