module EnsembleGeneMap

function ensembleToHugo()
    ensembletohugo = Dict{}()

    for line in readlines("./data/mart_export.txt")
        linefields = split(chomp(line), '\t')
        if linefields[2] != ""
            ensembletohugo[linefields[1]] = linefields[2]
        end
    end

    return ensembletohugo
end

end
