module Inference

#include("probability.jl")
#include("evidence.jl")
#include("nlmodel.jl")
#include("results.jl")
#include("keyoutputs.jl")
#include("expression.jl")
#include("pi.jl")
include("idMap.jl")

function run(config, verbose)
    println(config);
    if haskey(config, "id-map")
        if verbose
            println("Reading in idMap")
        end
        IDMap = IdMap.get(config["id-map"])
        println("got map")
        exit()
    else
        prinlnt("Need to specify 'id-map' in config file")
        exit(1)
    end
    println(IDMap)
    exit()
#    expression = Expression.get(IDMap, config["tissueType"])


#    for pathway in config["pathways"]
#        runPathway(config, pathway, expression, IdMap, verbose)
#    end
end

#function runPathway(config, pathway, expression, IDMap, verbose)
#    pinodes = Pi.readFile(pathway)
#
#
#    observations = Observations.get(config["evidence"], pinodes, IDMap)
#    keyoutputs = Keyoutputs.getNodes(config["keyoutputs"])
#    keyoutputs = ()
#    nodesampleresults = Dict()
#    allNodeResults = Dict()
#    for sample in observations["columns"]
#        if sample == "gene"
#            continue
#        end
#
#        if verbose
#            println("Running $sample")
#        end
#
#        samplenodestate = observations["samplenodestate"]
#        nodestate = samplenodestate[sample]
#
#        (sampleresults, x, x_bar) = NLmodel.run(pinodes,
#                                                nodestate,
#                                                keyoutputs,
#                                                lowerbound,
#                                                upperbound,
#                                                expression,
#                                                verbose)
#
#        for nodeName in keys(sampleresults)
#            if length(keys(nodesampleresults)) == 0 || haskey(nodesampleresults, nodeName) == false
#                nodesampleresults[nodeName] = Dict()
#            end
#            nodesampleresults[nodeName][sample] = sampleresults[nodeName]
#        end
#
#        if allOutputsFile != ""
#           for nodeName in keys(x)
#               if length(keys(allNodeResults)) == 0 || haskey(allNodeResults, nodeName) == false
#                   allNodeResults[nodeName] = Dict()
#               end
#               allNodeResults[nodeName][sample] = [x[nodeName], x_bar[nodeName]]
#           end
#        end
#    end
#    Results.createcsv(nodesampleresults,
#                      observations["columns"],
#                      resultsfile)
#
#    if allOutputsFile != ""
#         Results.outputAllResults(allNodeResults,
#                                  observations["columns"],
#                                  allOutputsFile)
#    end
#end
#
#function analyzeResults(resultsfile, expectedfile, downregulatedcutoff, upregulatedcutoff, pgmlab, verbose)
#    expected_data = Results.getExpected(expectedfile)
#    expected = expected_data["samplenodestate"]
#    expected_counts = expected_data["counts"]
#
#    results_data = Results.getResults(resultsfile, downregulatedcutoff, upregulatedcutoff, pgmlab)
#    results = results_data["samplenodestate"]
#    results_counts = results_data["counts"]
#    probability = 1;
#    correct_counts = Dict("1" => 0, "2" => 0, "3" => 0)
#    errors = Dict()
#    for patientname in keys(expected)
#        expected_patient_nodes = expected[patientname]
#        results_patient_nodes = results[patientname]
#        for nodename in keys(expected_patient_nodes)
#            expectedvalue = expected_patient_nodes[nodename]
#            resultvalue = results_patient_nodes[nodename]
#            if expectedvalue == resultvalue
#                correct_counts[expectedvalue] = correct_counts[expectedvalue] + 1
#            else
#                if haskey(errors, patientname) == false
#                    errors[patientname] = Dict()
#                end
#                errors[patientname][nodename] = [expectedvalue, resultvalue];
#            end
#        end
#    end
#    println("expected")
#    println(expected_counts)
#    println("results")
#    println(results_counts)
#    println("correct counts")
#    println(correct_counts)
#
#    total = 0
#    for (state, value) in expected_counts
#        total += value
#    end
#
#    total_correct = 0
#    for (state, value) in correct_counts
#        total_correct += value
#    end
#
#    percent_correct = total_correct/total * 100
#    println("percent correct: $percent_correct")
#
#    total_prob = Float64(1.0)
#    for (state, value) in correct_counts
#        new_prob = Probability.KorMoreSuccess(expected_counts[state], value, expected_counts[state]/total)
#        total_prob *= new_prob
#    end
#
#    println("Probability: $total_prob\n")
#
#    println("Patient\tNode\tExpected\tActual\n")
#    for patientname in keys(errors)
#        patient = errors[patientname]
#        for genename in keys(patient)
#            values = errors[patientname][genename]
#            println("$patientname\t$genename\t", values[1], "\t", values[2])
#        end
#    end
#end

end
