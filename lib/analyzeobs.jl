module analyzeObs

function run(pi_nodes)
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
