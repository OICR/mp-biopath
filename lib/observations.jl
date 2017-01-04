module Observations

include("dbidnamemapping.jl")

function testIdxs(fname, nodes)
    nodesList = collect(keys(nodes))

    measuredIdxs = Array{Integer}()

    if ismatch(r"test", fname)
        measuredIdxs = indexin(["b", "a"], nodesList)
    elseif ismatch(r"Signaling_by_ERBB2", fname)
        measuredIdxs = indexin(["1963571", "p-10Y-ERBB3-1_[plasma_membrane]_54424"], nodesList)
    elseif ismatch(r"DNA_Double-Strand_Break_Repair", fname)
        measuredIdxs = indexin(["RAD52_[nucleoplasm]_62640",
                       #"PALB2_[nucleoplasm]_241569",
                       "BRCA2_[nucleoplasm]_50952"
                    ], nodesList)
    elseif ismatch(r"PIP3_activates_AKT_signaling", fname)
        #test AKT upregulated
        #measuredidxs = indexin(["AKT1_[cytosol]_58253", "AKT2_[cytosol]_49860", "AKT3_[cytosol]_415917"], nodesList)

        #test PIK3CA upgregulated test
        #measuredidxs = indexin(["PIK3CA_[cytosol]_61074"], nodesList)

        #test3 PTEN downregulated
        measuredIdxs = indexin(["PTEN_mRNA_[cytosol]_2318745"], nodesList)
    elseif ismatch(r"DNA_Damage_Reversal", fname)
        measuredIdxs = indexin(["5657646"], nodesList)
    end

    return measuredIdxs
end

function copynumberIdxs(fname, pinodes, dbidfile)
    genetonodes = DbIdNameMapping.geneToRootNodes(pinodes, dbidfile)

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
                        #divide by two because it is copy number data
                       samplenodestate[column][node] = parse(Float64,lineparts[i])/ 2
                    end
                end
            end
        end
    end

    return Dict("columns" => headerparts, "samplenodestate" => samplenodestate)
end

end
