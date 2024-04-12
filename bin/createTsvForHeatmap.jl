#!/usr/bin/env julia

using ArgParse

include("../src/createTsvForHeatmap.jl")

function parse_commandline()
    s = ArgParseSettings("This scrpt creates a TSV file in the format required for creating a heatmap",
                         version = "1.0.1",
                         add_version = true)

    @add_arg_table s begin
        "--verbose", "-v"
            help = "This will cause output to be printed to standard out."
            action = :store_true
        "results-folder"
            help = "This is the directory that the inference results for a run are located."
            required = true
        "pathway-list-file"
            help = "This is the pathway list file provided by Reactome"
            required = true
        "key-outputs-file"
            help = "This file contains the list of keyoutputs for each pathway"
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

    CreateTsvForHeatmap.run(parsed_args["results-folder"],
                            parsed_args["pathway-list-file"],
                            parsed_args["key-outputs-file"],
                            parsed_args["verbose"])
end

main()
