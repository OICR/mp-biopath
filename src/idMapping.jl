module DbIdNameMapping

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
