module IdMap

function getIDmap(IDfile)
    columns = [:Database_Identifier,
               :Node_Name,
               :Node_Type,
               :Display_Name,
               :Reference_Entity_Name,
               :Reference_Entity_Identifier,
               :Instance_Class]

    idMap = Dict()
    open(IDfile, "r") do file
        header = readline(file)
        for line in eachline(file)
            lineParts = split(chomp(line), "\t")
            row = Dict(zip(columns, lineParts))
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

            refName = row[:Reference_Entity_Name]
            if refName != "N/A"
                if haskey(idMap, refName) == false
                    idMap[refName] = []
                end
                push!(idMap[refName], row)
            end

            refID = row[:Reference_Entity_Identifier]
            if refID != "N/A"
                if haskey(idMap, refID) == false
                    idMap[refID] = []
                end
                push!(idMap[refID], row)
            end
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
