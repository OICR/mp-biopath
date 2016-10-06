measuredIdxs = Array{Integer}()

if ismatch(r"test", ARGS[1])
    measuredIdxs = indexin(["b", "a"], nodesList)
elseif ismatch(r"Signaling_by_ERBB2", ARGS[1])
    measuredIdxs = indexin(["1963571", "p-10Y-ERBB3-1_[plasma_membrane]_54424"], nodesList)
elseif ismatch(r"DNA_Double-Strand_Break_Repair", ARGS[1])
    measuredIdxs = indexin(["RAD52_[nucleoplasm]_62640",
                   #"PALB2_[nucleoplasm]_241569",
                   "BRCA2_[nucleoplasm]_50952"
                ], nodesList)
elseif ismatch(r"PIP3_activates_AKT_signaling", ARGS[1])
    #test AKT upregulated
    #measuredidxs = indexin(["AKT1_[cytosol]_58253", "AKT2_[cytosol]_49860", "AKT3_[cytosol]_415917"], nodesList)

    #test PIK3CA upgregulated test
    #measuredidxs = indexin(["PIK3CA_[cytosol]_61074"], nodesList)

    #test3 PTEN downregulated
    measuredIdxs = indexin(["PTEN_mRNA_[cytosol]_2318745"], nodesList)
elseif ismatch(r"DNA_Damage_Reversal", ARGS[1])
    measuredIdxs = indexin(["5657646"], nodesList)
end

