module Keyoutputs

using CSV
using DataFrames
using DataStructures

function getNodes(keyoutputsFilename)

    columns = OrderedDict("parent_pathway_id" => String,
                 "pathway_id" => String,
                 "pathway_name" => String,
                 "output_reactome_dbid" => String,
                 "mapping_id" => String,
                 "description" => String,
                 "specific" => String,
                 "direct_measurement" => String,
                 "dm_method" => String,
                 "indirect_measurement" => String,
                 "im_method" => String)

    df = CSV.read(keyoutputsFilename,
                  delim = "\t",
                  header = collect(keys(columns)),
                  types = collect(values(columns)),
                  datarow = 2)


    keyoutputs = Dict{AbstractString,Any}()
    for row in eachrow(df)
        mapping_id = row[:mapping_id]
        if (isnull(mapping_id) == false)
            keyoutputs[get(mapping_id)] = 1
        end
    end

    return Set(collect(keys(keyoutputs)))
end

end
