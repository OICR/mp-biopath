module Essential

include("observations.jl")

function ensembleToHugo()
    ensembletohugo = Dict{}()

    for line in readlines("./data/mart_export2.txt")
        linefields = split(chomp(line), '\t')
        if linefields[2] != ""
            ensembletohugo[linefields[1]] = linefields[2]
        end
    end

    return ensembletohugo
end

function essentialGenesHugo()
    ensembletohugo = ensembleToHugo()
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


function getNodes()
    essentialgenesHugo = essentialGenesHugo()
    genetonode = Observations.geneToNodes()

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
