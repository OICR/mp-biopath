module Inference

include("probability.jl")
include("evidence.jl")
include("nlmodel.jl")
include("results.jl")
include("expression.jl")
include("pi.jl")
include("idMap.jl")

function run(configFile, config, verbose)
    if haskey(config, "id-map")
        if verbose
            println("Reading in idMap")
        end
        IDMap = IdMap.get(config["id-map"])
    else
        prinln("Need to specify 'id-map' in config file")
        exit(1)
    end

    if haskey(config, "expression")
        if verbose
            println("Reading in expression")
        end

        if haskey(config["expression"], "file") == false || haskey(config["expression"], "tissue") == false
            println("Expresion section of config needs 'file' and 'tissue' to be specified")
            exit(1)
        end

        expression = Expression.getTissue(config["expression"]["file"],
                                          config["expression"]["tissue"])

    else
        println("Missing expression section of config")
        exit(1)
    end

    if haskey(config, "upperbound")
        upperbound = config["upperbound"]
    else
        println("Need to specify upperbound")
        exit(1)
    end

    if haskey(config, "lowerbound")
        lowerbound = config["lowerbound"]
    else
        println("Need to specify lowerbound")
        exit(1)
    end

    runID = haskey(config, "id")? config["id"]: Base.Random.uuid4()
    info("runID: $runID")

    if haskey(config, "results")
        resultsConfig = config["results"]
        if haskey(resultsConfig, "directory")

             directory = resultsConfig["directory"]
             runDir = endswith(directory, "/")? "$directory$runID": "$directory/$runID"
        else
            println("Need to specify directory in results section of config")
            exit(1)
        end
    else
        println("Results section required")
        exit(1)
    end

    mkpath(runDir)
    cp(configFile, "$runDir/config.yaml"; remove_destination=true)

    if haskey(config, "evidence")
        evidence = Evidence.getEvidence(config["evidence"], IDMap, runDir)
    else
        println("Evidence section of config is missing")
        exit(1)
    end


    for file in config["pathways"]
        if verbose
            println("Running pathway: $file")
        end

        m = match(r".*/(?<pathway>.*)\.tsv", file)
        pathwayName = m[:pathway]
        pathwayDir = "$runDir/$pathwayName"
        runPathway(file, expression, IdMap, evidence, lowerbound, upperbound, config["coin-options"], pathwayDir, verbose)
    end
end

function runPathway(file, expression, IDMap, evidence, lowerbound, upperbound, options, pathwayDir, verbose)
    pinodes = Pi.readFile(file)
    nodesampleresults = Dict()
    allNodeResults = Dict()
    for sample in keys(evidence)
        if verbose
            println("Running $sample")
        end

        nodestate = evidence[sample]

        (sampleresults, x, x_bar) = NLmodel.runModel(pinodes,
                                                    nodestate,
                                                    lowerbound,
                                                    upperbound,
                                                    expression,
                                                    options,
                                                    verbose)

        for nodeName in keys(sampleresults)
            if length(keys(nodesampleresults)) == 0 || haskey(nodesampleresults, nodeName) == false
                nodesampleresults[nodeName] = Dict()
            end
            nodesampleresults[nodeName][sample] = sampleresults[nodeName]
        end

        for nodeName in keys(x)
            if length(keys(allNodeResults)) == 0 || haskey(allNodeResults, nodeName) == false
                allNodeResults[nodeName] = Dict()
            end
            allNodeResults[nodeName][sample] = [x[nodeName], x_bar[nodeName]]
        end
    end

    mkpath(pathwayDir)
    resultsPath = "$pathwayDir/results.tsv"

    if verbose
        println("Outputing results: $resultsPath")
    end

    Results.createcsv(nodesampleresults,
                      keys(evidence),
                      resultsPath)

end

function analyzeResults(resultsfile, expectedfile, downregulatedcutoff, upregulatedcutoff, pgmlab, verbose)
    expected_data = Results.getExpected(expectedfile)
    expected = expected_data["samplenodestate"]
    expected_counts = expected_data["counts"]

    results_data = Results.getResults(resultsfile, downregulatedcutoff, upregulatedcutoff, pgmlab)
    results = results_data["samplenodestate"]
    results_counts = results_data["counts"]
    probability = 1;
    correct_counts = Dict("1" => 0, "2" => 0, "3" => 0)
    errors = Dict()
    for patientname in keys(expected)
        expected_patient_nodes = expected[patientname]
        results_patient_nodes = results[patientname]
        for nodename in keys(expected_patient_nodes)
            expectedvalue = expected_patient_nodes[nodename]
            resultvalue = results_patient_nodes[nodename]
            if expectedvalue == resultvalue
                correct_counts[expectedvalue] = correct_counts[expectedvalue] + 1
            else
                if haskey(errors, patientname) == false
                    errors[patientname] = Dict()
                end
                errors[patientname][nodename] = [expectedvalue, resultvalue];
            end
        end
    end
    println("expected")
    println(expected_counts)
    println("results")
    println(results_counts)
    println("correct counts")
    println(correct_counts)

    total = 0
    for (state, value) in expected_counts
        total += value
    end

    total_correct = 0
    for (state, value) in correct_counts
        total_correct += value
    end

    percent_correct = total_correct/total * 100
    println("percent correct: $percent_correct")

    total_prob = Float64(1.0)
    for (state, value) in correct_counts
        new_prob = Probability.KorMoreSuccess(expected_counts[state], value, expected_counts[state]/total)
        total_prob *= new_prob
    end

    println("Probability: $total_prob\n")

    println("Patient\tNode\tExpected\tActual\n")
    for patientname in keys(errors)
        patient = errors[patientname]
        for genename in keys(patient)
            values = errors[patientname][genename]
            println("$patientname\t$genename\t", values[1], "\t", values[2])
        end
    end
end

end
