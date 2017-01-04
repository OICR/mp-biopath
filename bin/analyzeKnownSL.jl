#!/usr/bin/env julia

using ArgParse

include("../lib/analyzeknownsllist.jl")

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
        "--verbose", "-v"
            help = "This will cause output to be printed to standard out."
            action = :store_true
        "pairwise-interaction-file"
            help = "This is the full path to the pairwise interaction file."
            required = true
        "db-id-file"
            help = "This is the full path to the pairwise interaction file."
            required = true
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

    AnalyzeKnownSLList.run(parsed_args["pairwise-interaction-file"],
                           parsed_args["lowerbound"],
                           parsed_args["upperbound"],
                           parsed_args["downregulated-cutoff"],
                           parsed_args["upregulated-cutoff"],
                           parsed_args["db-id-file"],
                           parsed_args["verbose"])
end

main()
