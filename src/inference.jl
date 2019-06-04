module Inference

using DataFrames

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
            println("Reading in idMap")
        end
        idMap = IdMap.getIDmap(config["id-map"])
    else
        println("ERROR: Need to specify 'id-map' in config file")
        exit(1)
    end

    tissue = Dict()
    expression = DataFrame()
    if haskey(config, "expression")
        if verbose
            println("Reading in expression")
        end

        if haskey(config, "expression")
            if haskey(config["expression"], "file") == false || haskey(config["expression"], "tissue") == false
                println("ERROR: Expresion section of config needs 'file' and 'tissue' to be specified")
                exit(1)
            else
                expression = Expression.getExpression(config["expression"]["file"])
                tissueParam = config["expression"]["tissue"]
                if haskey(tissueParam, "mapping-file")
                    tissue["map"] = SampleTissueMap.getTissueMap(tissueParam["mapping-file"])
                elseif haskey(tissueParam, "name")
                    tissue["name"] = tissueParam["name"]
                else
                    println("ERROR: Under the expression -> tissue section of the config you need to specify either the key mapping-file or name")
                    exit(1)
                end
            end
        end
    end

    runID = haskey(config, "id") ? config["id"] : Base.Random.uuid4()
    println("runID: $runID")

    if haskey(config, "outputs")
        resultsConfig = config["outputs"]
        if haskey(resultsConfig, "directory")
             dir = resultsConfig["directory"]
             resultsDir = endswith(dir, "/") ? "$dir$runID" : "$dir/$runID"
             pathwayResultsDir = "$resultsDir/pathways"
        else
            error("Need to specify directory in results section of config")
        end
    else
        error("Results section in config required")
    end

    mkpath(resultsDir)
    mkpath(pathwayResultsDir)
    cp(configFile, "$resultsDir/config.yaml", force=true)

    if haskey(config, "evidence")
        evidence = Evidence.getEvidence(config["evidence"], idMap, resultsDir)
    else
        error("Evidence section of config is missing")
    end

    if @isdefined(expression)
        expressionColumns = names(expression)
        for sample in expressionColumns
            if haskey(tissue, "name")
                tissueType = Symbol(tissue["name"])
            elseif haskey(tissue, "map")
                if haskey(tissue["map"], Symbol(sample))
                    tissueType = Symbol(tissue["map"][Symbol(sample)])
                else
                    println("Tissue for sample $sample not found")
                    exit(1)
                end
            end

            if @isdefined(tissueType) && findfirst(isequal(tissueType),expressionColumns) == 0
                println("ERROR tissue type $tissueType not found in expression file")
                exit(1)
            end
        end
    end

    for file in config["pathways"]
        if verbose
            println("Running pathway: $file")
        end

        m = match(r".*/(?<pathway>.*)\.tsv", file)
        pathwayName = m[:pathway]
        pathwayDir = "$pathwayResultsDir/$pathwayName"
        runPathway(file, expression, idMap, evidence, tissue, pathwayDir, verbose)
    end
end

function isValueNormal(k,v)::Bool
  if v == 1.0
      return true
  else
      return false
  end
end

function runPathway(file, expression, IDMap, evidence, tissue, pathwayDir, verbose)
    pinodes = Pi.readFile(file)
    nodeSampleResults = Dict()
    basePathwayTissueResults = Dict()
    for sample in keys(evidence)
        if haskey(tissue, "name") || haskey(tissue, "map")
            if haskey(tissue, "name")
                tissueType = Symbol(tissue["name"])
            elseif haskey(tissue, "map")
                if haskey(tissue["map"], Symbol(sample))
                    tissueType = Symbol(tissue["map"][Symbol(sample)])
                else
                    println("Tissue for sample $sample not found")
                    exit(1)
                end
            end

            expressionMap = Expression.getTissue(expression, tissueType)
        else
             tissueType = "unspecified"
             expressionMap = Dict()
        end

        nodeState = evidence[sample]
        println(nodeState)
        if haskey(basePathwayTissueResults, tissueType) == false
            controlNodeState = Dict()
            for (k, v) in nodeState
               if v == 1.0
                  controlNodeState[k] = v
               end
            end
            println(controlNodeState)
            (basePathwayResults, x, x_bar) = NLmodel.runModel(pinodes,
                                                              controlNodeState,
                                                              expressionMap,
                                                              verbose)
            basePathwayTissueResults[tissueType] = basePathwayResults
        end

        if verbose
            println("Running $sample")
        end

        (sampleResults, x, x_bar) = NLmodel.runModel(pinodes,
                                                     nodeState,
                                                     expressionMap,
                                                     verbose)

        for nodeName in keys(sampleResults)
            if length(keys(nodeSampleResults)) == 0 || haskey(nodeSampleResults, nodeName) == false
                nodeSampleResults[nodeName] = Dict()
            end

            calculatedValue = sampleResults[nodeName] - (basePathwayTissueResults[tissueType][nodeName] - 1)
            nodeSampleResults[nodeName][sample] = (calculatedValue < 0.01) ? 0.01 : calculatedValue
        end
    end

    mkpath(pathwayDir)
    resultsPath = "$pathwayDir/results.tsv"

    if verbose
        println("Outputing results: $resultsPath")
    end

    Results.createcsv(nodeSampleResults,
                      keys(evidence),
                      resultsPath)
end

end
