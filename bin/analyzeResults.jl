#!/usr/bin/env julia

using ArgParse

include("../lib/inference.jl")

function parse_commandline()
    s = ArgParseSettings("This interface is for analyzing the results files produced by running inference with respect to a expected results file.",
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
            default = 1.01
        "--upperbound", "-u"
            arg_type = Int
            default = 10
        "--lowerbound", "-l"
            arg_type = Float64
            default = 0.001
        "--pgmlab"
            help = "This is for makeing it so the numbers for the states are converted to 1,2,3 for down nc and up"
            action = :store_true
        "--verbose", "-v"
            help = "This will cause output to be printed to standard out."
            action = :store_true
        "results-file"
            help = "This file will be generated and contain the resulting inferred data."
            required = true
        "expected-file"
            help = "This is the corresponding expected state file to the observation file."
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
                             parsed_args["pgmlab"],
                             parsed_args["verbose"])
end

main()
