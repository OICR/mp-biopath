module Observations

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

function geneToNode()
    genetonode = Dict()

    for line in readlines("data/db_id_to_name_mapping.txt")
        lineparts = split(chomp(line), "\t")
        if match(r"Reference", lineparts[3]) != nothing
            if in(genetonode, lineparts[2])
                push(genetonode[lineparts[2]], lineparts[1])
            else
                genetonode[lineparts[2]] = [lineparts[1]]
            end
        end
    end

    return genetonode
end

function copynumberIdxs(fname)

    genetonode = geneToNode()

    data = readlines(fname)

    header = shift!(data)
    headerparts = split(chomp(header), "\t")

    samplenodestate = Dict()

    i = 0
    for column in headerparts
        i = i + 1
        if i != 1
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
                for node in genetonode[gene]
                    samplenodestate[column][node] = parse(Float64,lineparts[i])
                end
            end
        end

    end



    return Dict("columns" => headerparts, "samplenodestate" => samplenodestate)
end

end