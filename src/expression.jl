module Expression

using CSV
using DataFrames
using DataStructures

include("dbidnamemapping.jl")

function hugoGeneExpression(tissueType)
    hugoGeneExpressionValues = Dict()

    expressionFilename = "./data/expression-data.tsv"

    columns = OrderedDict("Gene ID" => String,
                   "Gene Name" => String,
                   "adipose tissue" => Float16,
                   "adrenal gland" => Float16,
                   "bone marrow" => Float16,
                   "cerebral cortex" => Float16,
                   "colon" => Float16,
                   "duodenum" => Float16,
                   "endometrium" => Float16,
                   "esophagus" => Float16,
                   "fallopian tube" => Float16,
                   "gall bladder" => Float16,
                   "heart" => Float16,
                   "kidney" => Float16,
                   "liver" => Float16,
                   "lung" => Float16,
                   "lymph node" => Float16,
                   "ovary" => Float16,
                   "pancreas" => Float16,
                   "placenta" => Float16,
                   "prostate gland" => Float16,
                   "rectum" => Float16,
                   "saliva-secreting gland" => Float16,
                   "skeletal muscle tissue" => Float16,
                   "small intestine" => Float16,
                   "smooth muscle tissue" => Float16,
                   "spleen" => Float16,
                   "stomach" => Float16,
                   "testis" => Float16,
                   "thyroid gland" => Float16,
                   "tonsil" => Float16,
                   "urinary bladder" => Float16,
                   "vermiform appendix" => Float16,
                   "zone of skin" => Float16)

    df = CSV.read(expressionFilename,
                  delim = "\t",
                  header = collect(keys(columns)),
                  types = collect(values(columns)),
                  datarow = 5)

    for row in eachrow(df)
        hugoGeneExpressionValues[row[Symbol("Gene Name")]] = row[Symbol(tissueType)]

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
