#!/usr/bin/env julia

using ArgParse

include("../lib/pi.jl")
include("../lib/observations.jl")
include("../lib/nlmodel.jl")
include("../lib/keyoutputs.jl")
include("../lib/essential.jl")

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

   # essentialgenes = Essential.getNodes()
    keyoutputs = Keyoutputs.getNodes()

    nodes = Pi.readFile(parsed_args["pairwise-interaction-file"])

    measuredIdxs = Array{Integer}()
    if parsed_args["find-si"]
        println("to do")
    else
        if parsed_args["observation-file"] != nothing
            observations = Observations.copynumberIdxs(parsed_args["observation-file"])
            for sample in observations["columns"]
                if sample == "Gene"
                    continue
                end
                if parsed_args["verbose"]
                    println("Running $sample")
                end
                samplenodestate = observations["samplenodestate"]
                nodestate = samplenodestate[sample]
                NLmodel.run(nodes,
                            nodestate,
                            parsed_args["lowerbound"],
                            parsed_args["upperbound"],
                            parsed_args["downregulated-cutoff"],
                            parsed_args["upregulated-cutoff"],
                            parsed_args["verbose"])

            end
        else
            #in future need to make it say a warning that you need to provide an observation file if --find-si not specified
            measuredIdxs = Observations.testIdxs(parsed_args["pairwise-interaction-file"], nodes)
            if parsed_args["verbose"]
                 println("Measured indexes: $measuredIdxs")
            end
        end
    end
    exit()
end

main()
