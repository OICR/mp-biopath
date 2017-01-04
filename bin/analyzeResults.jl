#!/usr/bin/env julia

using ArgParse

include("../lib/inference.jl")

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
        "results-file"
            help = "This is the full path to the results file. (required if providing observation file)"
            required = true
        "expected-file"
            help = "this is the corresponding expected state file to the observation file. Used when analyzing results"
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

    Inference.analyzeResults(parsed_args["results-file"],
                             parsed_args["expected-file"],
                             parsed_args["downregulated-cutoff"],
                             parsed_args["upregulated-cutoff"],
                             parsed_args["verbose"])
end

main()
