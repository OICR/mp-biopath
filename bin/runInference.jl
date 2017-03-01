#!/usr/bin/env julia

using ArgParse

include("../lib/inference.jl")

function parse_commandline()
    s = ArgParseSettings("This script is for infering values for keyoutputs for pathway from observation data.",
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
        "--copynumber"
        help = "This flag will make it so the values of the nodes are divided by two"
        action = :store_true
        "--verbose", "-v"
        help = "This will cause output to be printed to standard out."
        action = :store_true
        "pairwise-interaction-file"
        help = "This is the full path to the pairwise interaction file."
        required = true
        "observation-file"
        help = "This is the full path to the observation file."
        required = true
        "results-file"
        help = "This is the full path to the results file generated by this script"
        required = true
        "db-id-file"
        help = "This is the full path to the db_id_to_name_mapping file."
        required = true
        "key-outputs-file"
        help = "This is the full path of the keyoutputs file."
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

    Inference.run(parsed_args["pairwise-interaction-file"],
                  parsed_args["observation-file"],
                  parsed_args["results-file"],
                  parsed_args["key-outputs-file"],
                  parsed_args["db-id-file"],
                  parsed_args["lowerbound"],
                  parsed_args["upperbound"],
                  parsed_args["downregulated-cutoff"],
                  parsed_args["upregulated-cutoff"],
                  parsed_args["copynumber"],
                  parsed_args["verbose"])
end

main()
