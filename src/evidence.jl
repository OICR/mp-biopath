module Evidence

using Nullables
using DataFrames
using CSV

function getEvidence(evidence, idMapping, runDir)
    evidenceMapFile = "$runDir/evidenceIDMapping.tsv"

    sampleNodeValue = Dict()
    if haskey(evidence, "dna")
        evidenceDNA = evidence["dna"]
        if haskey(evidenceDNA, "genomic")
            genomicFile = evidenceDNA["genomic"]
            genomicEvidence = getGenomic(genomicFile, idMapping)
        end
        sampleNodeValue = genomicEvidence["sampleNodeValue"]
        outputToEvidenceMap("genomic", evidenceMapFile, genomicEvidence["geneNodesMap"])
    end

    return sampleNodeValue
end

function outputToEvidenceMap(evidenceType, file, geneNodesMap)
    outfile = open(file, "w")

    for gene in keys(geneNodesMap)
        for node in geneNodesMap[gene]
            id = node[:Database_Identifier]
            nodeName = node[:Node_Name]
            write(outfile, "$evidenceType\t$gene\t$id\t$nodeName\n")
        end
    end

    close(outfile)
end

function getGenomic(file, idMap)
    headers = split(strip(readline(file)),'\t')
    ncols = length(headers)
    datatypes = vcat(String, [Float64 for i=2:ncols])

    df = CSV.read(file, delim="\t", types=datatypes, weakrefstrings=false)

    geneNodeMap = Dict()
    sampleNodeValue = Dict()
    geneNodesMap = Dict()
    for row in eachrow(df)
        gene = get(row[Symbol("gene")])
        if haskey(idMap, gene) == true
            nodes = idMap[gene]
            geneNodesMap[gene] = nodes;
            geneValue = Dict()
            first = true
            for entry in row
                if first == false
                    sample = entry[1]
                    if isnull(entry[2]) == false
                        value = get(entry[2])
                        if value != -999
                            if haskey(sampleNodeValue, sample) == false
                                sampleNodeValue[sample] = Dict()
                            end
                            for node in nodes
                                sampleNodeValue[sample][node[:Database_Identifier]] = (value == 0)? 0.01: value
                            end
                        end
                    end
                end
                first = false
            end
        end
    end
    return Dict("sampleNodeValue" => sampleNodeValue,
                "geneNodesMap" => geneNodesMap)
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
