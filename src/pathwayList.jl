module PathwayList

using CSV
using DataFrames

function getPathwayList(file)
    df = CSV.read(file,
                  delim="\t",
                  datarow=1,
                  header=["pathway_id",
                          "pathway_name",
                          "short_name",
                          "colour"],
                  types=[String,
                         String,
                         String,
                         String])
    pathwayList = Dict()
    for row in eachrow(df)
        if !isnothing(row[:pathway_name])
            pathwayName = row[:pathway_name]
            pathwayList[pathwayName] = row
        end
    end
  
    return pathwayList
end

function getPathwayNameToIdMap(file)
    df = CSV.read(file,
                  delim="\t",
                  datarow=1,
                  header=["pathway_id",
                          "pathway_name",
                          "short_name",
                          "colour"],
                  types=[String,
                         String,
                         String,
                         String])
    pathwayList = Dict()
    for row in eachrow(df)
        if !isnothing(row[:pathway_name])
            pathwayName = row[:pathway_name]
            pathwayList[pathwayName] = row[:pathway_id]
        end
    end
  
    return pathwayList
end



end
