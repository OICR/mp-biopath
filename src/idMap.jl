module IdMap


using DataFrames
using Nullables
using CSV

  #                null="N/A",
function get(IDfile)
    df = CSV.read(IDfile,
                  delim="\t",
                  datarow=2,
                  quotechar="\\",
                  nullable=false,
                  header=["Database_Identifier",
                          "Node_Name",
                          "Display_Name",
                          "Reference_Entity_Name",
                          "Reference_Entity_Identifier",
                          "Instance_Class"],
                  types=[String,
                         String,
                         String,
                         String,
                         String,
                         String])

    idMap = Dict()
    for row in eachrow(df)
        dbID = row[:Database_Identifier]
        if haskey(idMap, dbID) == false
            idMap[dbID] = []
        end
        push!(idMap[dbID], row)

        nodeName = row[:Node_Name]
        if haskey(idMap, nodeName) == false
            idMap[nodeName] = []
        end
        push!(idMap[nodeName], row)

        refID = row[:Reference_Entity_Identifier]
        if refID != "N/A"
            if haskey(idMap, refID) == false
                idMap[refID] = []
            end
            push!(idMap[refID], row)
        end
    end

    return idMap
end

function idToNodes(idMap, id, entityType="generic")
   if entityType == "generic"
       println("need to finish this")
   end

   return 1
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
