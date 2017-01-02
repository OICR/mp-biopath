#!/usr/bin/env julia

using ArgParse

include("../lib/pi.jl")
include("../lib/analyzeobs.jl")
include("../lib/findsi.jl")
include("../lib/analyzeknownsilist.jl")

function parse_commandline()
    s = ArgParseSettings("This program infers the value of nodes in Reactome pathways from observation data.",
                         version = "0.0.1",
                         add_version = true)

    @add_arg_table s begin
        "--downregulated-cutoff"
            help = "This determines at which level the node is determined to be down regulated."
            arg_type = Float64
            default = 0.99
        "--upregulated-cutoff"
            help = "This determines at which level the node is determined to be upregulated."
            arg_type = Float64
            default = 1.02
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
        "results-file"
        help = "This is the full path to the results file. (required if providing observation file)"
            required = false
        "expected-file"
            help = "this is the corresponding expected state file to the observation file. Used when analyzing results"
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
        FindSI.run(pinodes,
                   parsed_args["lowerbound"],
                   parsed_args["upperbound"],
                   parsed_args["downregulated-cutoff"],
                   parsed_args["upregulated-cutoff"],
                   parsed_args["pairwise-interaction-file"],
                   parsed_args["verbose"])
    elseif parsed_args["analyze-known-si-list"]
        AnalyzeKnownSiList.run(pinodes,
                               parsed_args["lowerbound"],
                               parsed_args["upperbound"],
                               parsed_args["downregulated-cutoff"],
                               parsed_args["upregulated-cutoff"],
                               parsed_args["pairwise-interaction-file"],
                               parsed_args["verbose"])
    else
        if parsed_args["expected-file"] != nothing
            AnalyzeObs.inspect(parsed_args["observation-file"],
                               parsed_args["expected-file"],
                               parsed_args["downregulated-cutoff"],
                               parsed_args["upregulated-cutoff"],
                               parsed_args["verbose"])
            println("after inspect")
        else
            AnalyzeObs.run(pinodes,
                           parsed_args["observation-file"],
                           parsed_args["results-file"],
                           parsed_args["key-outputs"],
                           parsed_args["lowerbound"],
                           parsed_args["upperbound"],
                           parsed_args["downregulated-cutoff"],
                           parsed_args["upregulated-cutoff"],
                           parsed_args["verbose"])
        end
    end
end

main()
