module SampleTissueMap

using CSV
using DataFrames

function getTissueMap(file)
    df = CSV.read(file, delim="\t", weakrefstrings=false, nullable=false)
    sampleTissueMap = Dict()
    for row in eachrow(df)
        sampleTissueMap[Symbol(row[Symbol("sample")])] = row[Symbol("tissue_type")]
    end

    return sampleTissueMap
end

end
