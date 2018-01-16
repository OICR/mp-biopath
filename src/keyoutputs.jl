module KeyOutputs

using CSV
using DataFrames

function getKeyoutputs(file)
    df = CSV.read(file,
                  delim="\t",
                  datarow=2,
                  quotechar="\\",
                  header=["pathway_id",
                          "pathway_name",
                          "output_reactome_dbid",
                          "mapping_id",
                          "description",
                          "specific",
                          "direct_measurement",
                          "dm_method",
                          "indirect_measurement",
                          "im_methodDatabase_Identifier"],
                  types=[String,
                         String,
                         String,
                         String,
                         String,
                         String,
                         String,
                         String,
                         String,
                         String])

    keyoutputMap = Dict()
    for row in eachrow(df)
        pathwayID = get(row[Symbol("pathway_id")])
        if haskey(keyoutputMap, pathwayID) == false
            keyoutputMap[pathwayID] = []
        end
        push!(keyoutputMap[pathwayID], get(row[:output_reactome_dbid]))
    end
  
    return keyoutputMap
end

end
