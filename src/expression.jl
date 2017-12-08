module Expression

using CSV
using DataFrames
using DataTables
using DataStructures

function geneExpressionMap(filename, dataRow, tissueType)
    # Where the Gene Name is the HUGO name
    # The data is expected to be Float16 data
    df = CSV.read(filename,
                  delim = "\t",
                  datarow = dataRow)

    geneExpressionValues = Dict()
    for row in eachrow(df)
        geneExpressionMap[row[Symbol("Gene Name")]] = row[Symbol(tissueType)]
    end

    return geneExpressionMap
end

end
