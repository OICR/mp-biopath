module Essential

include("observations.jl")

function getNodes()
    ensembleToHugo = Dict{}()

    for line in readlines("./data/mart_export2.txt")
        linefields = split(chomp(line), '\t')
        if linefields[2] != ""
            ensembleToHugo[linefields[1]] = linefields[2]
        end
    end

    essentialgenesHugo = ASCIIString[]

    for line in readlines("./data/essential_9606_all_gene_ids.txt")
        ensembleGene = chomp(line)
        if haskey(ensembleToHugo, ensembleGene)
            push!(essentialgenesHugo, ensembleToHugo[ensembleGene])
        end
    end



    genetonode = Observations.geneToNode()

    essentialnodes = ASCIIString[]

    for hugogene in essentialgenesHugo
        if haskey(genetonode, hugogene)
            for node in genetonode[hugogene]
                push!(essentialnodes, node)
            end
        end
    end

    return essentialnodes
end

function getGenes(nodes)
    return collect(intersect(Set(nodes), Set(getNodes())))
end

end
