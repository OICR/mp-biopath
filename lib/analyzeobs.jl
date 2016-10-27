module AnalyzeObs

using ExcelReaders

include("observations.jl")
include("nlmodel.jl")
include("results.jl")

function run(pinodes, observationfile, keyoutputsfile, lowerbound, upperbound, downregulatedcutoff, upregulatedcutoff, verbose)
    observations = Observations.copynumberIdxs(observationfile, pinodes)
    keyoutputs = keyoutputsfile? Keyoutputs.getNodes(): Set()

    nodesampleresults = Dict()
    for sample in observations["columns"]
        if sample == "Gene"
            continue
        end

        if verbose
            println("Running $sample")
        end

        samplenodestate = observations["samplenodestate"]
        nodestate = samplenodestate[sample]

        sampleresults = NLmodel.run(pinodes,
                                    nodestate,
                                    keyoutputs,
                                    lowerbound,
                                    upperbound,
                                    downregulatedcutoff,
                                    upregulatedcutoff,
                                    verbose)

        for nodeName in keys(sampleresults)
            if length(keys(nodesampleresults)) == 0 || haskey(nodesampleresults, nodeName) == false
                nodesampleresults[nodeName] = Dict()
            end
            nodesampleresults[nodeName][sample] = sampleresults[nodeName]
        end
    end
    Results.createcsv(nodesampleresults,
                      observations["columns"],
                      observationfile)
end


function inspect(observationfile, expectedfile, downregulatedcutoff, upregulatedcutoff, verbose)

    expected = Results.getExpected(expectedfile)

    results = Results.getResults(observationfile, downregulatedcutoff, upregulatedcutoff)


    println(results)

end

end
