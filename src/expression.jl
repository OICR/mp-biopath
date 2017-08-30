module Expression

include("dbidnamemapping.jl")

function hugoGeneExpression(tissueType)
    hugoGeneExpressionValues = Dict()

    column_index = -1
    line_count = 0
    for line in readlines("./data/expression-data.tsv")
        if ismatch(r"^#", line) == false
            line_count += 1
            fields = chomp(line)
            line_parts = split(fields, "\t")

            if line_count == 1
                column_index = findfirst(x -> x == tissueType, line_parts)
                if column_index == 0
                    println("tissue type not found")
                    exit()
                end
            else
                if column_index == -1
                    exit()
                else
                    hugoGene = convert(String, line_parts[2])
                    expression_value = line_parts[column_index]
                    hugoGeneExpressionValues[hugoGene] = expression_value
                end
            end
        end
    end

    return hugoGeneExpressionValues
end


function getNodes(dbidfile, tissueType)
    hugoGeneExpressionValues = hugoGeneExpression(tissueType)
    genetonode = DbIdNameMapping.geneToNodes(dbidfile)

    expressionNodeValues = Dict()

    for (hugogene, value) in hugoGeneExpressionValues
        if haskey(genetonode, hugogene)
            for node in genetonode[hugogene]
                expressionNodeValues[node] = value
            end
        end
    end

    return expressionNodeValues
end

function get(dbidfile, tissueType)
    return getNodes(dbidfile, tissueType)
end

end
