module Observations

function get(evidence, pinodes, IDMapping)
    sampleNodeState = Dict()

    if evidence["dna"]
        getGenomic(evidence["dna"], IDMapping)
    end


end

function getDnaEvidence(evidence, idMapping)
    divideby = copynumber? 2:1
    data = readlines(fname)

    header = shift!(data)
    headerparts = split(chomp(header), "\t")

    samplenodestate = Dict()

    i = 0
    for column in headerparts
        i = i + 1
        if i > 1
            samplenodestate[column] = Dict()
        end
    end

    gene = ""
    for line in data
        lineparts = split(chomp(line), "\t")
        i = 0
        for column in headerparts
            i = i + 1
            if i == 1
                gene = lineparts[1]
            else
                if haskey(genetonodes, gene)
                    for node in Array(genetonodes[gene])
                        samplenodestate[column][node] = oneNormal?
                                                            parse(Float64, lineparts[i]) :
                                                            copynumber?
                                                               parse(Float64,lineparts[i]) / 2 :
                                                               parse(Float64,lineparts[i]) - 1
                    end
                end
            end
        end
    end

    return Dict("columns" => headerparts, "samplenodestate" => samplenodestate)
end

end
