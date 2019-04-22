module Inference

include("probability.jl")
include("evidence.jl")
include("nlmodel.jl")
include("results.jl")
include("expression.jl")
include("pi.jl")
include("idMap.jl")
include("sampleTissueMap.jl")

function run(configFile, config, verbose)

    if haskey(config, "id-map")
        if verbose
            info("Reading in idMap")
        end
        IDMap = IdMap.getIDmap(config["id-map"])
    else
        info("ERROR: Need to specify 'id-map' in config file")
        exit(1)
    end

    expression = Dict()
    if haskey(config, "expression")
        if verbose
            info("Reading in expression")
        end

        if haskey(config["expression"], "file") == false || haskey(config["expression"], "tissue") == false
            info("ERROR: Expresion section of config needs 'file' and 'tissue' to be specified")
            exit(1)
        end

        expression = Expression.getExpression(config["expression"]["file"])
        tissueParam = config["expression"]["tissue"]

        tissue = Dict()
        if haskey(tissueParam, "mapping-file")
            tissue["map"] = SampleTissueMap.getTissueMap(tissueParam["mapping-file"])
        elseif haskey(tissueParam, "name")
            tissue["name"] = tissueParam["name"]
        else
            info("ERROR: Under the expression -> tissue section of the config you need to speficy either the key mapping-file or name")
            exit(1)
        end
    end

    if haskey(config, "upperbound")
        upperbound = config["upperbound"]
    else
        info("ERROR: Need to specify upperbound in config")
        error(1)
    end

    if haskey(config, "lowerbound")
        lowerbound = config["lowerbound"]
    else
        warn("Need to specify lowerbound in config")
        exit(1)
    end

    runID = haskey(config, "id")? config["id"]: Base.Random.uuid4()
    info("runID: $runID")

    if haskey(config, "results")
        resultsConfig = config["results"]
        if haskey(resultsConfig, "directory")
             directory = resultsConfig["directory"]
             resultsDir = endswith(directory, "/")? "$directory$runID": "$directory/$runID"
             pathwayResultsDir = "$resultsDir/pathways"
        else
            error("Need to specify directory in results section of config")
        end
    else
        error("Results section in config required")
    end

    mkpath(resultsDir)
    mkpath(pathwayResultsDir)
    cp(configFile, "$resultsDir/config.yaml"; remove_destination=true)

    if haskey(config, "evidence")
        evidence = Evidence.getEvidence(config["evidence"], IDMap, resultsDir)
    else
        error("Evidence section of config is missing")
    end

    for sample in keys(evidence)
        if haskey(tissue, "name")
            tissueType = Symbol(tissue["name"])
        else
            if haskey(tissue["map"], Symbol(sample))
                tissueType = Symbol(tissue["map"][Symbol(sample)])
            else
                info("Tissue for sample $sample not found")
                exit(1)
            end
         end

         expressionColumns = names(expression)
         if findfirst(expressionColumns, tissueType) == 0
	     info("ERROR tissue type $tissueType not found in expression file")
             exit(1)
         end
    end

    @sync begin
        for file in config["pathways"]
            if verbose
                println("Running pathway: $file")
            end

            m = match(r".*/(?<pathway>.*)\.tsv", file)
            pathwayName = m[:pathway]
            pathwayDir = "$pathwayResultsDir/$pathwayName"
            @async runPathway(file, expression, IdMap, evidence, lowerbound, upperbound, config["coin-options"], tissue, pathwayDir, verbose)
        end
    end
end

function runPathway(file, expression, IDMap, evidence, lowerbound, upperbound, options, tissue, pathwayDir, verbose)
    pinodes = Pi.readFile(file)
    nodeSampleResults = Dict()
    basePathwayTissueResults = Dict()

    for sample in keys(evidence)
        if haskey(tissue, "name")
            tissueType = Symbol(tissue["name"])
        else
            if haskey(tissue["map"], Symbol(sample))
                tissueType = Symbol(tissue["map"][Symbol(sample)])
            else
                info("Tissue for sample $sample not found")
		exit(1)
            end
        end
        expressionMap = tissueType == Expression.getTissue(expression, tissueType)

        if haskey(basePathwayTissueResults, tissueType) == false
            (basePathwayResults, x, x_bar) = NLmodel.runModel(pinodes,
                                                              Dict(),
                                                              lowerbound,
                                                              upperbound,
                                                              expressionMap,
                                                              options,
                                                              verbose)
            basePathwayTissueResults[tissueType] = basePathwayResults
        end

        if verbose
            info("Running $sample")
        end

        nodeState = evidence[sample]

        (sampleResults, x, x_bar) = NLmodel.runModel(pinodes,
                                                     nodeState,
                                                     lowerbound,
                                                     upperbound,
                                                     expressionMap,
                                                     options,
                                                     verbose)

        for nodeName in keys(sampleResults)
            if length(keys(nodeSampleResults)) == 0 || haskey(nodeSampleResults, nodeName) == false
                nodeSampleResults[nodeName] = Dict()
            end
	    calculatedValue = sampleResults[nodeName] - (basePathwayTissueResults[tissueType][nodeName] - 1)
            nodeSampleResults[nodeName][sample] = (calculatedValue < 0) ? 0 : calculatedValue
        end
    end

    mkpath(pathwayDir)
    resultsPath = "$pathwayDir/results.tsv"

    if verbose
        info("Outputing results: $resultsPath")
    end

    Results.createcsv(nodeSampleResults,
                      keys(evidence),
                      resultsPath)
end

