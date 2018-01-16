module PathwayList

using CSV
using DataFrames

function getPathwayList(file)
    df = CSV.read(file,
                  delim="\t",
                  datarow=1,
                  quotechar="\\",
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
        if !isnull(row[Symbol("pathway_name")])
            pathwayName = get(row[Symbol("pathway_name")])
            pathwayList[pathwayName] = row
        end
    end
  
    return pathwayList
end

end
