module Expression

using CSV
using DataFrames

function getTissue(file, tissue)
    # Where the Gene Name is the HUGO name
    df = CSV.read(file, delim="\t", weakrefstrings=false)
    tissueValues = Dict()
    for row in eachrow(df)
        value = row[Symbol(tissue)]
        if isnull(value) == false
            tissueValues[get(row[Symbol("Gene Name")])] = get(value)
        end
    end

    return tissueValues
end

end
