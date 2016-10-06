#!/usr/bin/env julia

using ArgParse

include("../lib/pi.jl")
include("../lib/model.jl")

function parse_commandline()
    s = ArgParseSettings("This program infers the value of nodes in Reactome pathways from observation data",
                         version = "0.0.1",
                         add_version = true)

    @add_arg_table s begin
        "--downregulated-cutoff"
            help = "This determines at which level the node is determined to be down regulated"
            arg_type = Float64
            default = 0.9
        "--upregulated-cutoff"
            help = "This determines at which level the node is determined to be upregulated"
            arg_type = Float64
            default = 1.1
        "--upperbound", "-u"
            arg_type = Int
            default = 10
        "--lowerbound", "-l"
            arg_type = Float64
            default = 0.001
        "--verbose", "-v"
            help = "This will cause output to be printed to standard out"
            action = :store_true
        "pairwise-interaction-file"
            help = "This is the full path to the pairwise interaction file"
            required = true
        "observation-file"
            help = "This is the full path to the observation file."
            required = false
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    #interactions = get_interactions()
    #measured_nodes = get_measured_nodes()
    #model = run_model()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
end

main()
