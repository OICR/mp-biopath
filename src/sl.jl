module SL

include("dbidnamemapping.jl")
include("observations.jl")

function getNodes(dbidfile)
    genetonodes = DbIdNameMapping.geneToNodes(dbidfile)

    genenodes = getGenes()
    for gene in genenodes
        geneasymbolnodes = haskey(genetonodes, gene["GeneASymbol"])? genetonodes[gene["GeneASymbol"]]: []
        genebsymbolnodes = haskey(genetonodes, gene["GeneBSymbol"])? genetonodes[gene["GeneBSymbol"]]: []
        gene["GeneANodes"] = geneasymbolnodes
        gene["GeneBNodes"] = genebsymbolnodes
    end

    return genenodes
end

function getGenes()
    genes = []
    header = true
    columns = Array{Dict{String,Any}}[]
    for line in readlines("data/sl_human")
        if header
            columns = split(chomp(line), "\t")
            header = false
        else
            lineparts = split(chomp(line), "\t")
            nodeone = String(lineparts[1])
            nodetwo = String(lineparts[3])
            push!(genes, Dict{String,Any}(zip(columns, lineparts)))
        end
    end

    return genes
end

end
