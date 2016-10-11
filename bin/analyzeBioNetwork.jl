#!/usr/bin/env julia

using ArgParse

include("../lib/pi.jl")
include("../lib/observations.jl")
include("../lib/nlmodel.jl")
include("../lib/keyoutputs.jl")
include("../lib/essential.jl")
include("../lib/results.jl")

function parse_commandline()
    s = ArgParseSettings("This program infers the value of nodes in Reactome pathways from observation data.",
                         version = "0.0.1",
                         add_version = true)

    @add_arg_table s begin
        "--downregulated-cutoff"
            help = "This determines at which level the node is determined to be down regulated."
            arg_type = Float64
            default = 0.5
        "--upregulated-cutoff"
            help = "This determines at which level the node is determined to be upregulated."
            arg_type = Float64
            default = 1.5
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

    nodes = Pi.readFile(parsed_args["pairwise-interaction-file"])

    if parsed_args["find-si"]
        essentialgenes = Essential.getGenes(collect(keys(nodes)))

        quasiessentialnodestates = Dict()

        for i in [0,2]
            for node in keys(nodes)
                if contains(node, "PSEUDONODE")
                    continue
                end
                if in(node, Set(essentialgenes)) == false || i == 2
                    sampleresults = NLmodel.run(nodes,
                                                Dict(node => i),
                                                essentialgenes,
                                                parsed_args["lowerbound"],
                                                parsed_args["upperbound"],
                                                parsed_args["downregulated-cutoff"],
                                                parsed_args["upregulated-cutoff"],
                                                parsed_args["verbose"])

                    for noderesult in values(sampleresults)
                        state = NLmodel.valueToState(noderesult,
                                                     parsed_args["downregulated-cutoff"],
                                                     parsed_args["upregulated-cutoff"])
                        if state == "Down Regulated"
                            quasiessentialnodestates[node] = i
                        end
                    end
                end
            end
        end

        pairwisefile = parsed_args["pairwise-interaction-file"]
        sifilename = join([pairwisefile, "si"], ".")
        sioutfile = open(sifilename, "w")
        write(sioutfile, "count\tnodeone\tnodeone_value\tnode_two\tnoed_two_value\teffected_node\teffected_node_value\n")
        flush(sioutfile)

        count = 0
        for i in [0,2]
            for j in [0,2]
                for nodeone in keys(nodes)
                    for nodetwo in keys(nodes)
                        if contains(nodeone, "PSEUDONODE") || contains(nodetwo, "PSEUDONODE")
                            continue
                        end

                        if (nodeone == nodetwo && i == j) || in(nodeone, Set(essentialgenes)) || in(nodetwo, Set(essentialgenes))
                            continue
                        end

                        if (haskey(quasiessentialnodestates, nodeone) && (quasiessentialnodestates[nodeone] == i))
                            continue
                        end

                        if (haskey(quasiessentialnodestates, nodetwo) && (quasiessentialnodestates[nodetwo] == j))
                            continue
                        end

                        sampleresults = NLmodel.run(nodes,
                                                    Dict(nodeone => i,
                                                         nodetwo => j),
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
                                count = count + 1
                                write(sioutfile, "$count\t$nodeone\t$i\t$nodetwo\t$j\t$resultnode\t$value\n")
                                flush(sioutfile)
                            end
                        end
                    end
                end
            end
        end
        close(sioutfile)
    else
        if parsed_args["observation-file"] != nothing
            observations = Observations.copynumberIdxs(parsed_args["observation-file"])
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

                sampleresults = NLmodel.run(nodes,
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
