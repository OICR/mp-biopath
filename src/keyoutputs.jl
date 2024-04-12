module KeyOutputs

using CSV
using DataFrames

function getKeyoutputs(file)
    df = CSV.read(file;
                  delim='\t',
                  quotechar='\\',
                  header=["pathway_id",
                          "pathway_name",
                          "node_id",
                          "node_name"],
                  types=[String,
                         String,
                         String,
                         String])
    keyoutputList = Dict()
    for row in eachrow(df)
        if !isnothing(row[:pathway_name])
            pathwayName = row[:pathway_name]
            if haskey(keyoutputList, pathwayName) == false
               keyoutputList[pathwayName] = Dict()
            end
            keyoutputList[pathwayName][row[:node_id]] = row[:node_name]
        end
    end
  
    return keyoutputList
end

end
