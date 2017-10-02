module IdMap

using DataFrames
using CSV

function get(dbidfile)
    df = CSV.read(fname,
                  delim="\t",
                  nullable=false,
                  header=["Database_Identifier",
                          "Node_Name",
                          "Display_Name",
                          "Reference_Entity_Name",
                          "Reference_Entity_Identifier",
                          "Instance_Class",
                          "Type"],
                  types=[int,
                         String,
                         String,
                         String,
                         String,
                         String,
                         String])

    idMapping = Dict()    
    for row in eachrow(df)
       nodeType = row[:Type]

       if haskey(idMapping, nodeType)
           idMapping[nodeType] = Dict()
       end

       entityID = row[:referenceEntityID]
       if haskey(idMapping[nodeType], entityID)
          idMapping[nodeType][entityID] = []
       end
       push!(idMapping[nodeType][entityID], row)

       entityName = row[:referenceEntityName]
       if haskey(idMapping[nodeType], entityName)
          idMapping[nodeType][entityName] = []
       end
       push!(idMapping[nodeType][entityName], row)
    end

    return idMapping
end


function allGeneReferenceProduct(dbidfile)
    genes = []
    header = true
    for line in readlines(dbidfile)
        if header == false
            lineparts = split(chomp(line), "\t")
            node = String(lineparts[2])
            if contains(lineparts[6], "Reference") == true
                push!(genes, node)
            end
        end
        header = false
    end

    return genes
end

function geneToNodes(dbidfile)
    genetonodes = Dict()
    header = true
    for line in readlines(dbidfile)
        if header == false
            lineparts = split(chomp(line), "\t")
            node = String(lineparts[2])
            gene = String(lineparts[3])
            #if match(r"Reference", lineparts[6]) != nothing
            if haskey(genetonodes, gene)
                nodes = genetonodes[gene]
                push!(nodes, node)
            else
                genetonodes[gene] = [node]
            end
        end
        header = false
    end

    return genetonodes
end

function nodeToGene(dbidfile)
   genetonodes = geneToNodes(dbidfile)

   nodetogene = Dict()
   for (gene, nodes) in genetonodes
       for node in nodes
           nodetogene[node] = gene
       end
   end

   return nodetogene
end

function geneToRootNodes(pinodes, dbidfile)
    genetonodes = Dict()
    header = true
    for line in readlines(dbidfile)
        if header == false
            lineparts = split(chomp(line), "\t")
            node = String(lineparts[2])
            gene = String(lineparts[3])
            if haskey(pinodes, node) && pinodes[node].relation == "ROOT"
    #        if match(r"Reference", lineparts[6]) != nothing && haskey(pinodes, node) && pinodes[node].relation == "ROOT"
                if haskey(genetonodes, gene)
                    nodes = genetonodes[gene]
                    push!(nodes, node)
                else
                    genetonodes[gene] = [node]
                end
            end
        end
        header = false
    end

    return genetonodes
end

end
