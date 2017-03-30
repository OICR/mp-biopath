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
    columns = Array{Dict{ASCIIString,Any}}[]
    for line in readlines("data/sl_human")
        if header
            columns = split(chomp(line), "\t")
            header = false
        else
            lineparts = split(chomp(line), "\t")
            nodeone = ASCIIString(lineparts[1])
            nodetwo = ASCIIString(lineparts[3])
            push!(genes, Dict{ASCIIString,Any}(zip(columns, lineparts)))
        end
    end

    return genes
end

end
