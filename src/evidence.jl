module Evidence

using Nullables
using DataFrames
using CSV

function getEvidence(evidence, idMapping)
    sampleNodeValue = Dict()
    if haskey(evidence, "dna")
        evidenceDNA = evidence["dna"]
        if haskey(evidenceDNA, "genomic")
            genomicFile = evidenceDNA["genomic"]
            genomicEvidence = getGenomic(genomicFile, idMapping)
        end
        sampleNodeValue = genomicEvidence["sampleNodeValue"]
    end

    return sampleNodeValue
end

function getGenomic(file, idMap)
    df = CSV.read(file, delim="\t", weakrefstrings=false)

    geneNodeMap = Dict()
    sampleNodeValue = Dict()
    for sample in eachrow(df)
        gene = get(sample[Symbol("gene")])
        nodes = idMap[gene]
        geneNodeMap[gene] = nodes

        geneValue = Dict()
        first = true
        for entry in sample
            if first == false
                sample = entry[1]
                if isnull(entry[2]) == false
                    if haskey(sampleNodeValue, sample) == false
                        sampleNodeValue[sample] = Dict()
                    end
                    value = get(entry[2])
                    for node in nodes
                        sampleNodeValue[sample][node[:Node_Name]] = value - 1
                    end
                end
            end
            first = false
        end
    end

    return Dict("sampleNodeValue" => sampleNodeValue,
                "geneNodeMap" => geneNodeMap)
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
