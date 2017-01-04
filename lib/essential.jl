module Essential

include("ensemblegenemap.jl")
include("dbidnamemapping.jl")

function essentialGenesHugo()
    ensembletohugo = EnsembleGeneMap.ensembleToHugo()
    essentialgenesHugo = ASCIIString[]

    for line in readlines("./data/essential_9606_all_gene_status.txt")
        fields = chomp(line)
        line_parts = split(fields)
        ensembleGene = line_parts[1]
        if haskey(ensembletohugo, ensembleGene)
            push!(essentialgenesHugo, ensembletohugo[ensembleGene])
        end
    end

    return essentialgenesHugo
end


function getNodes(dbidfile)
    essentialgenesHugo = essentialGenesHugo()
    genetonode = DbIdNameMapping.geneToNodes(dbidfile)

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

function getGenes(nodes, dbidfile)
    return collect(intersect(Set(nodes), Set(getNodes(dbidfile))))
end

end
