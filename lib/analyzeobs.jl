module AnalyzeObs

using ExcelReaders

include("observations.jl")
include("nlmodel.jl")
include("results.jl")
include("keyoutputs.jl")

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

    expected_data = Results.getExpected(expectedfile)


    expected = expected_data["samplenodestate"]
    expected_counts = expected_data["counts"]

    results_data = Results.getResults(observationfile, downregulatedcutoff, upregulatedcutoff)
    results = results_data["samplenodestate"]
    results_counts = results_data["counts"]
    probability = 1;
    correct_counts = Dict("1" => 0, "2" => 0, "3" => 0)
    for patientname in keys(expected)
        expected_patient_nodes = expected[patientname]
        results_patient_nodes = results[patientname]
        for nodename in keys(expected_patient_nodes)
            expectedvalue = expected_patient_nodes[nodename]
            resultvalue = results_patient_nodes[nodename]
            if expectedvalue == resultvalue
                correct_counts[expectedvalue] = correct_counts[expectedvalue] + 1
            end
        end
    end
    println(correct_counts)
    println("done")
end

end
