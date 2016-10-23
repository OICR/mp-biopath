#!/usr/bin/env julia

using ArgParse

include("../lib/pi.jl")
include("../lib/observations.jl")
include("../lib/nlmodel.jl")
include("../lib/keyoutputs.jl")
include("../lib/essential.jl")
include("../lib/results.jl")
include("../lib/sl.jl")


function parse_commandline()
    s = ArgParseSettings("This program infers the value of nodes in Reactome pathways from observation data.",
                         version = "0.0.1",
                         add_version = true)

    @add_arg_table s begin
        "--downregulated-cutoff"
            help = "This determines at which level the node is determined to be down regulated."
            arg_type = Float64
            default = 0.9
        "--upregulated-cutoff"
            help = "This determines at which level the node is determined to be upregulated."
            arg_type = Float64
            default = 1.1
        "--upperbound", "-u"
            arg_type = Int
            default = 10
        "--lowerbound", "-l"
            arg_type = Float64
            default = 0.001
        "--find-si"
            help = "When this option is set do not provide any observations. This will systemattically analyze the network and find all synthetically lethal pairs."
            action = :store_true
        "--key-outputs"
            help = "If this is set it will prepare output based on ./data/keyoutputs.tsv"
            action = :store_true
        "--essential-genes"
            help = "If this is specified a report will be made with regards to essential genes"
        "--analyze-known-si-list"
            help = "This will go through sl_human and produce a file containing resulting values for the essential genes"
            action = :store_true
        "--verbose", "-v"
            help = "This will cause output to be printed to standard out."
            action = :store_true
        "pairwise-interaction-file"
            help = "This is the full path to the pairwise interaction file."
            required = true
        "observation-file"
            help = "This is the full path to the observation file."
            required = false
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    if parsed_args["verbose"]
        println("Parsed args:")
        for (arg,val) in parsed_args
            println("  $arg  =>  $val")
        end
    end

    pinodes = Pi.readFile(parsed_args["pairwise-interaction-file"])

    if parsed_args["find-si"]
        essentialgenes = Essential.getGenes(collect(keys(pinodes)))
        allGenes = Observations.allGeneReferenceProduct()
        allGenesSet = Set(allGenes)

        downregulatednodestates = []
        for node in keys(pinodes)
            if in(node, allGenesSet) == false
                continue
            end
            if contains(node, "PSEUDONODE") || in(node, Set(essentialgenes))
                continue
            end

            for i in [0,2]
                if in(node, Set(essentialgenes)) == false || i == 2
                    sampleresults = NLmodel.run(pinodes,
                                                Dict(node => i),
                                                essentialgenes,
                                                parsed_args["lowerbound"],
                                                parsed_args["upperbound"],
                                                parsed_args["downregulated-cutoff"],
                                                parsed_args["upregulated-cutoff"],
                                                parsed_args["verbose"])

                    essentialnodes = []
                    for resultnode in keys(sampleresults)
                        value = sampleresults[resultnode]
                        state = NLmodel.valueToState(value,
                                                     parsed_args["downregulated-cutoff"],
                                                     parsed_args["upregulated-cutoff"])

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

        pairwisefile = parsed_args["pairwise-interaction-file"]
        sifilename = join([pairwisefile, "si"], ".")
        sloutfile = open(sifilename, "w")
        write(sloutfile, "count\tnodeone\tnodeone_value\tnodetwo\tnodetwo_value\teffected_node\teffected_node_value\n")
        flush(sloutfile)

        count = 0
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
                candidates = candidates + 1

                sampleresults = NLmodel.run(pinodes,
                                            Dict(nodeone["name"] => nodeone["state"],
                                                 nodetwo["name"] => nodetwo["state"]),
                                                 essentialgenes,
                                                 parsed_args["lowerbound"],
                                                 parsed_args["upperbound"],
                                                 parsed_args["downregulated-cutoff"],
                                                 parsed_args["upregulated-cutoff"],
                                                 parsed_args["verbose"])

                for resultnode in keys(sampleresults)
                   value = sampleresults[resultnode]
                   state = NLmodel.valueToState(value,
                                                parsed_args["downregulated-cutoff"],
                                                parsed_args["upregulated-cutoff"])

                   if state == "Down Regulated"
                       for essentialnodeone in nodeone["essentialNodes"]
                            if essentialnodeone["name"] != resultnode
                                continue
                            end
                            for essentialnodetwo in nodetwo["essentialNodes"]
                                if essentialnodetwo["name"] != resultnode
                                    continue
                                end

                                count = count + 1
                                write(sloutfile, "$count\t")
                                write(sloutfile, nodeone["name"])
                                write(sloutfile, "\t")
                                write(sloutfile, essentialnodeone["value"])
                                write(sloutfile, "\t")
                                write(sloutfile, nodetwo["name"])
                                write(sloutfile, "\t")
                                write(sloutfile, essentialnodetwo["value"])
                                write(sloutfile, "\t")
                                write(sloutfile, "$resultnode\t$value\n")
                                flush(sloutfile)

                                break
                            end
                            break
                        end
                    end
                end
            end
        end
        println("number of candidate SL $candidates\n")
        close(sloutfile)

    elseif parsed_args["analyze-known-si-list"]
        slnodes = SL.getNodes()
        pinodesSet = Set(keys(pinodes))

        essentialgenes = Essential.getGenes(collect(keys(pinodes)))
        nodetogene     = Observations.nodeToGene()

        slessential = Dict{ASCIIString,Any}[]
        slpinodes = ASCIIString[]
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
        allGenes = Observations.allGeneReferenceProduct()
        allGenesSet = Set(allGenes)
        for node in slpinodes
            sampleresults = NLmodel.run(pinodes,
                                        Dict(node => 0),
                                        essentialgenes,
                                        parsed_args["lowerbound"],
                                        parsed_args["upperbound"],
                                        parsed_args["downregulated-cutoff"],
                                        parsed_args["upregulated-cutoff"],
                                        parsed_args["verbose"])

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

        pairwisefile = parsed_args["pairwise-interaction-file"]
        slfilename = join([pairwisefile, "si.analysis"], ".")
        sloutfile = open(slfilename, "w")
        headercolumns = ["Count"]
        for column in keys(slessential[1])
            push!(headercolumns, column)
        end
        header = join(headercolumns, "\t")
        write(sloutfile, string(header, "\n"))
        flush(sloutfile)

        for slnodespair in slpinodesPairs
            sampleresults = NLmodel.run(pinodes,
                                        Dict(slnodespair[1] => 0,
                                             slnodespair[2] => 0),
                                        essentialgenes,
                                        parsed_args["lowerbound"],
                                        parsed_args["upperbound"],
                                        parsed_args["downregulated-cutoff"],
                                        parsed_args["upregulated-cutoff"],
                                        parsed_args["verbose"])

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
            columnvalues = ASCIIString[]
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
    else
        if parsed_args["observation-file"] != nothing
            observations = Observations.copynumberIdxs(parsed_args["observation-file"], pinodes)
            keyoutputs = parsed_args["key-outputs"]? Keyoutputs.getNodes(): Set()

            nodesampleresults = Dict()
            for sample in observations["columns"]
                if sample == "Gene"
                    continue
                end

                if parsed_args["verbose"]
                    println("Running $sample")
                end

                samplenodestate = observations["samplenodestate"]
                nodestate = samplenodestate[sample]

                sampleresults = NLmodel.run(pinodes,
                                            nodestate,
                                            keyoutputs,
                                            parsed_args["lowerbound"],
                                            parsed_args["upperbound"],
                                            parsed_args["downregulated-cutoff"],
                                            parsed_args["upregulated-cutoff"],
                                            parsed_args["verbose"])

                for nodeName in keys(sampleresults)
                    if length(keys(nodesampleresults)) == 0 || haskey(nodesampleresults, nodeName) == false
                        nodesampleresults[nodeName] = Dict()
                    end
                    nodesampleresults[nodeName][sample] = sampleresults[nodeName]
                end
            end
            Results.createcsv(nodesampleresults,
                              observations["columns"],
                              parsed_args["observation-file"])
        end
    end
end

main()
