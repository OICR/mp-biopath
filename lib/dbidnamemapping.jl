module DbIdNameMapping

function allGeneReferenceProduct()
    genes = []
    for line in readlines("data/db_id_to_name_mapping.txt")
        lineparts = split(chomp(line), "\t")
        node = ASCIIString(lineparts[1])
        if contains(lineparts[3], "Reference") == true
            push!(genes, node)
        end
    end

    return genes
end

function geneToNodes()
    genetonodes = Dict()
    for line in readlines("data/db_id_to_name_mapping.txt")
        lineparts = split(chomp(line), "\t")
        node = ASCIIString(lineparts[1])
        gene = ASCIIString(lineparts[2])
        if match(r"Reference", lineparts[3]) != nothing
            if haskey(genetonodes, gene)
                nodes = genetonodes[gene]
                push!(nodes, node)
            else
                genetonodes[gene] = [node]
            end
        end
    end

    return genetonodes
end

function nodeToGene()
   genetonodes = geneToNodes()

   nodetogene = Dict()
   for (gene, nodes) in genetonodes
       for node in nodes
           nodetogene[node] = gene
       end
   end

   return nodetogene
end

function geneToRootNodes(pinodes)
    genetonodes = Dict()
    for line in readlines("data/db_id_to_name_mapping.txt")
        lineparts = split(chomp(line), "\t")
        node = ASCIIString(lineparts[1])
        gene = ASCIIString(lineparts[2])
        if match(r"Reference", lineparts[3]) != nothing && haskey(pinodes, node) && pinodes[node].relation == "ROOT"
            if haskey(genetonodes, gene)
                nodes = genetonodes[gene]
                push!(nodes, node)
            else
                genetonodes[gene] = [node]
            end
        end
    end

    return genetonodes
end

end
