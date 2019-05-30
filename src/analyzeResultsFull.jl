module AnalyzeResultsFull

using DataFrames
include("valuetostate.jl")
include("probability.jl")
include("evidence.jl")
include("results.jl")
include("idMap.jl")
include("pathwayList.jl")
include("pi.jl")

function addTotalsToConfusionMatrix(confusion_matrix)
    confusion_matrix[1,4] = confusion_matrix[1,1] + confusion_matrix[1,2] + confusion_matrix[1,3] #first row total
    confusion_matrix[2,4] = confusion_matrix[2,1] + confusion_matrix[2,2] + confusion_matrix[2,3]
    confusion_matrix[3,4] = confusion_matrix[3,1] + confusion_matrix[3,2] + confusion_matrix[3,3]
    

    confusion_matrix[4,1] = confusion_matrix[1,1] + confusion_matrix[2,1] + confusion_matrix[3,1] #first column total
    confusion_matrix[4,2] = confusion_matrix[1,2] + confusion_matrix[2,2] + confusion_matrix[3,2]
    confusion_matrix[4,3] = confusion_matrix[1,3] + confusion_matrix[2,3] + confusion_matrix[3,3]

    confusion_matrix[4,4] = confusion_matrix[1,4] + confusion_matrix[2,4] + confusion_matrix[3,4] #total of totals
  
    return confusion_matrix
end

function writeMatrix(matrix, filepath)
    open(filepath, "w") do f
         for row in eachrow(matrix)
             write(f, join(row, "\t"))
             write(f, "\n")
         end
    end
end

function addConfusionMatrixToDataFrame(df, cutoff, confusion_matrix)
   TP = confusion_matrix[1,1]
   FP = confusion_matrix[1,2]
   FN = confusion_matrix[2,1]
   TN = confusion_matrix[2,2]
   total = TN + FP + FN + TP
   P = TP + FN 
   N = TN + FP

   TPR = 0
   if P != 0
       TPR = TP / P
   end

   TNR = 0
   if N != 0
       TNR = TN / N
   end

   PPV = 0
   if (TP + FP) != 0
       PPV = TP / (TP + FP)
   end 

   NPV = 0
   if (TN + FN) != 0
       NPV = TN / (TN + FN)
   end

   FNR = 1 - TPR
   FPR = 1 - TNR
   FDR = 1 - PPV
   FOR = 1 - NPV
   ACC = (TP + TN) / total

   Fone = 0
   if (2 * TP + FP + FN) != 0
       Fone = (2 * TP) / (2 * TP + FP + FN)
   end
   row = [cutoff, TP, FP, FN, TN, total, P, N, TPR, TNR, PPV, NPV, FNR, FPR, FDR, FOR, ACC, Fone]
   #println(row)
   #println(dFloat64#print(length(row))
   push!(df, row)
end

function initializeConfusionMatricies()
   return DataFrame(Cutoff = Float64[],
                    TP = Int16[],
		    FP = Int16[],
		    FN = Int16[],
		    TN = Int16[],
		    Total = Int16[],
		    P = Int16[],
		    N = Int16[],
		    TPR = Float64[],
		    TNR = Float64[],
		    PPV = Float64[],
		    NPV = Float64[],
		    FNR = Float64[],
		    FPR = Float64[],
		    FDR = Float64[],
		    FOR = Float64[],
		    ACC = Float64[],
		    F1 = Float64[])
end

function analyzeResultsFull(results_folder, expected_results_folder, pathway_list_file, pathways_folder, db_id_name_mapping_file, analysis_results_folder, experimental_results_folder, verbose)
    if isfile(pathway_list_file)
       pathway_name_to_id_map = PathwayList.getPathwayNameToIdMap(pathway_list_file)
    else
       println("ERROR pathway list file does not exists: $pathway_list")
    end

    if isfile(db_id_name_mapping_file)
         id_map = IdMap.getIDmap(db_id_name_mapping_file)
    end

    expected_results = Dict()
    if isdir(expected_results_folder)
        for expected_results_filename in readdir(expected_results_folder)
            if endswith(expected_results_filename, "_expected_results.tsv")
                expected_results_filepath = joinpath(expected_results_folder, expected_results_filename)
                pathway_expected_results = Results.getExpected(expected_results_filepath)
                pathway_name = expected_results_filename[1:end - 21]
                pathway_id = pathway_name_to_id_map[pathway_name]
                expected_results[pathway_id] = pathway_expected_results
             end
        end
    else
        println("ERROR: Folder does not exist $expected_results_folder")
    end

    experimental_results = Dict()
    if isdir(experimental_results_folder)
        for experimental_results_filename in readdir(experimental_results_folder)
            if endswith(experimental_results_filename, "_experimental_results.tsv")
                experimental_results_filepath = joinpath(experimental_results_folder, experimental_results_filename)
                pathway_experimental_results = Results.getExpected(experimental_results_filepath)
                pathway_name = experimental_results_filename[1:end - 25]
                pathway_id = pathway_name_to_id_map[pathway_name]
                experimental_results[pathway_id] = pathway_experimental_results
             end
        end
    else
        println("ERROR: Folder does not exist $experimental_results_folder")
    end
    
    mp_biopath_results = Dict()
    if isdir(results_folder)
         main_tests_folder = joinpath(results_folder, "main-tests")
         if isdir(main_tests_folder)
             results_pathways_folder = joinpath(main_tests_folder, "pathways")
             if isdir(results_pathways_folder)
                 results_pathways_folders = readdir(results_pathways_folder)
                 for results_pathway_folder in results_pathways_folders
                       pathway_folder = joinpath(results_pathways_folder, results_pathway_folder)
                       pathway_id = pathway_name_to_id_map[results_pathway_folder]
                       pathway_results_filepath = joinpath(pathway_folder, "results.tsv")
                       if isfile(pathway_results_filepath)
                           mp_biopath_results[pathway_id] = Results.getResults(pathway_results_filepath)
                       else
                           println("ERROR: File does not exist: $pathway_results_filepath")
                           exit(4)
                       end
                 end

             else
                 println("ERROR Results pathway folder does not exits: $results_pathways_folder")
                 exit(6)
             end
         else
             println("ERROR Folder does not exist $main_tests_folder")
             exit(2)
         end
    else
         println("ERROR Folder Does not Exist: $results_folder")
         exit(1)
    end

    pathway_networks = Dict()
    if isdir(pathways_folder)
        for pi_file in readdir(pathways_folder)
             if endswith(pi_file, ".tsv")
                 pathway_name = pi_file[1:end - 4]
                 pathway_id = pathway_name_to_id_map[pathway_name]
                 pathway_networks[pathway_id] = Pi.readFile(joinpath(pathways_folder, pi_file))
             end
        end
    else
        println("ERROR: pathway folder does not exist: $pathways_folder")
        exit(7)
    end

    total_zeros_expected = 0
    total_ones_expected = 0
    total_twos_expected = 0
    for pathway_id in keys(expected_results)
        for scenario in keys(expected_results[pathway_id]["samplenodestate"])
            for keyoutput in keys(expected_results[pathway_id]["samplenodestate"][scenario])
                value = expected_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                if value == "0" 
                    total_zeros_expected += 1
                elseif value == "1"
                    total_ones_expected += 1
                elseif value == "2"
                    total_twos_expected += 1
                else
                    println(value)
                    println(typeof(value))
                    exit()
                end
            end
        end
    end

    total_tests = total_zeros_expected + total_ones_expected + total_twos_expected
    
    confusion_matrix = zeros(Int16, (4, 4))
    for pathway_id in keys(experimental_results)
        for scenario in keys(experimental_results[pathway_id]["samplenodestate"])
            for keyoutput in keys(experimental_results[pathway_id]["samplenodestate"][scenario])
                  expected_value = expected_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                  experimental_value = experimental_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                  if expected_value == "0" &&experimental_value == "0"
                      confusion_matrix[1,1] = confusion_matrix[1,1] + 1
                  elseif expected_value == "1" && experimental_value == "0"
                      confusion_matrix[2,1] = confusion_matrix[2,1] + 1
                  elseif expected_value == "2" && experimental_value == "0"
                      confusion_matrix[3,1] = confusion_matrix[3,1] + 1
                  elseif expected_value == "0" && experimental_value == "1"
                      confusion_matrix[1,2] = confusion_matrix[1,2] + 1
                  elseif expected_value == "1" && experimental_value == "1"
                      confusion_matrix[2,2] = confusion_matrix[2,2] + 1
                  elseif expected_value == "2" && experimental_value == "1"
                      confusion_matrix[3,2] = confusion_matrix[3,2] + 1
                  elseif expected_value == "0" && experimental_value == "2"
                      confusion_matrix[1,3] = confusion_matrix[1,3] + 1
                  elseif expected_value == "1" && experimental_value == "2"
                      confusion_matrix[2,3] = confusion_matrix[2,3] + 1
                  elseif expected_value == "2" && experimental_value == "2"
                      confusion_matrix[3,3] = confusion_matrix[3,3] + 1
                  end
             end
         end
    end
    
    row_totals = sum(confusion_matrix, dims=2)
    percent_zeros_expected = total_zeros_expected/row_totals[1] 
    percent_ones_expected = total_ones_expected/row_totals[2] 
    percent_twos_expected = total_twos_expected/row_totals[3] 
    
    biological_confusion_matrix = zeros(Float32, (4, 4))
    biological_confusion_matrix[1,1] = confusion_matrix[1,1] * percent_zeros_expected
    biological_confusion_matrix[1,2] = confusion_matrix[1,2] * percent_zeros_expected
    biological_confusion_matrix[1,3] = confusion_matrix[1,3] * percent_zeros_expected
    biological_confusion_matrix[2,1] = confusion_matrix[2,1] * percent_ones_expected
    biological_confusion_matrix[2,2] = confusion_matrix[2,2] * percent_ones_expected
    biological_confusion_matrix[2,3] = confusion_matrix[2,3] * percent_ones_expected
    biological_confusion_matrix[3,1] = confusion_matrix[3,1] * percent_twos_expected
    biological_confusion_matrix[3,2] = confusion_matrix[3,2] * percent_twos_expected
    biological_confusion_matrix[3,3] = confusion_matrix[3,3] * percent_twos_expected


    addTotalsToConfusionMatrix(confusion_matrix)
    addTotalsToConfusionMatrix(biological_confusion_matrix)

    writeMatrix(confusion_matrix, joinpath(analysis_results_folder, "expectedVsExperimentalConfusionMatrix.tsv"))
    writeMatrix(biological_confusion_matrix, joinpath(analysis_results_folder, "expectedVsExperimentalBiologicalConfusionMatrix.tsv"))

    ### determine optimal lowerbound cutoff based on F1 score
    lowerbound_confusion_matricies = initializeConfusionMatricies()
    
    for lowerbound in range(0,1, step=0.01)
        confusion_matrix = zeros(Int16, (2, 2))
        for pathway_id in keys(experimental_results)
            for scenario in keys(experimental_results[pathway_id]["samplenodestate"])
                for keyoutput in keys(experimental_results[pathway_id]["samplenodestate"][scenario])
                    expected_value = expected_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                    results_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
                    if expected_value == "2" || expected_value == "1"
                        if results_value > lowerbound
                            confusion_matrix[2,2] = confusion_matrix[2,2] + 1
                        else
                            confusion_matrix[2,1] = confusion_matrix[2,1] + 1
                        end
                    else 
                        if results_value > lowerbound
                            confusion_matrix[1,2] = confusion_matrix[1,2] + 1
                        else
                            confusion_matrix[1,1] = confusion_matrix[1,1] + 1
                        end
                    end
                end
            end
        end
        addConfusionMatrixToDataFrame(lowerbound_confusion_matricies, lowerbound, confusion_matrix)
    end
    println(lowerbound_confusion_matricies)

    upperbound_confusion_matricies = initializeConfusionMatricies()
    for upperbound in range(1,10, step=0.01)
        confusion_matrix = zeros(UInt16, (2, 2))
        for pathway_id in keys(experimental_results)
            for scenario in keys(experimental_results[pathway_id]["samplenodestate"])
                for keyoutput in keys(experimental_results[pathway_id]["samplenodestate"][scenario])
                    expected_value = expected_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                    results_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
                    if expected_value == "0" || expected_value == "1"
                        if results_value < upperbound
                            confusion_matrix[2,2] = confusion_matrix[2,2] + 1
                        else
                            confusion_matrix[2,1] = confusion_matrix[2,1] + 1
                        end
                    else 
                        if results_value < upperbound
                            confusion_matrix[1,2] = confusion_matrix[1,2] + 1
                        else
                            confusion_matrix[1,1] = confusion_matrix[1,1] + 1
                        end
                    end
                end
            end
        end
        addConfusionMatrixToDataFrame(upperbound_confusion_matricies, upperbound, confusion_matrix)
    end
    println(upperbound_confusion_matricies)
    exit()


    println(upperbound_confusion_matricies)
    exit()
    for index in 1:100
        println(index)
        lowerbound = 1 - index/100
        println(lowerbound)
    end


#    expected_data = Results.getExpected(expectedfile)
#    expected = expected_data["samplenodestate"]
#
#    expected_zero_got_zero = 0
#    expected_one_got_one = 0
#    expected_two_got_two = 0
#    expected_zero_got_one = 0
#    expected_zero_got_two = 0
#    expected_one_got_zero = 0
#    expected_one_got_two = 0
#    expected_two_got_one = 0
#    expected_two_got_zero = 0
#
#    results = Results.getResults(resultsfile)
#    probability = 1;
#    errors = Dict()
#    if full == true
#        println("scenario\tkeyoutput\tpredictiveValue")
#    end
#    for patientname in keys(expected)
#        expected_patient_nodes = expected[patientname]
#        results_patient_nodes = results[patientname]
#        for nodename in keys(expected_patient_nodes)
#            expectedState = expected_patient_nodes[nodename]
#            resultValue = results_patient_nodes[nodename]
#
#            resultState = ValueToState.getStateNumber(resultValue, downregulatedcutoff, upregulatedcutoff)
#            if expectedState != resultState
#                if haskey(errors, patientname) == false
#                    errors[patientname] = Dict()
#                end
#                errors[patientname][nodename] = [expectedState, resultValue]
#            end
#
#            if expectedState == "0" && resultState == "0"
#                if full == true
#                    println("$patientname\t$nodename\tTP")
#                else
#                    expected_zero_got_zero += 1
#                end
#            elseif expectedState == "1" && resultState == "1"
#                if full == true
#                    println("$patientname\t$nodename\tTN")
#                else
#                    expected_one_got_one += 1
#                end
#            elseif expectedState == "2" && resultState == "2"
#                if full == true
#                    println("$patientname\t$nodename\tTP")
#                else
#                    expected_two_got_two += 1
#                end
#            elseif expectedState == "0" && resultState == "1"
#                if full == true
#                    println("$patientname\t$nodename\tFN")
#                else
#                    expected_zero_got_one += 1
#                end
#            elseif expectedState == "0" && resultState == "2"
#                if full == true
#                    println("$patientname\t$nodename\tFP")
#                else
#                    expected_zero_got_two += 1
#                end
#            elseif expectedState == "1" && resultState == "0"
#                if full == true
#                    println("$patientname\t$nodename\tFP")
#                else
#                    expected_one_got_zero += 1
#                end
#            elseif expectedState == "1" && resultState == "2"
#                if full == true
#                    println("$patientname\t$nodename\tFP")
#                else
#                    expected_one_got_two += 1
#                end
#            elseif expectedState == "2" && resultState == "1"
#                if full == true
#                    println("$patientname\t$nodename\tFN")
#                else
#                    expected_two_got_one += 1
#                end
#            elseif expectedState == "2" && resultState == "0"
#                if full == true
#                    println("$patientname\t$nodename\tFP")
#                else
#                     expected_two_got_zero += 1
#                end
#            end
#        end
#    end
#
#    if full == true
#        return
#    end
#
#    TN = expected_one_got_one
#    println("\ntrue negative (TN): $TN")
#
#    TP = expected_two_got_two + expected_zero_got_zero
#    println("true positive (TP): $TP")
#
#    FP = expected_one_got_zero + expected_one_got_two + expected_two_got_zero + expected_zero_got_two
#    println("false positive (FP): $FP")
#
#    FN = expected_two_got_one + expected_zero_got_one
#    println("false negative (FN): $FN")
#
#    total = TN + FP + FN + TP
#    println("Total number of cases: $total")
#
#    P = TP + FN 
#    println("conditional positive (P): $P")    
#
#    N = TN + FP
#    println("conditional negative (N): $N")
#
#    TPR = 0
#    if P != 0
#        TPR = TP / P
#    end
#    println("sensitivity, recall, hit rate, or true positive rate (TPR): $TPR")
#
#    TNR = 0
#    if N != 0
#        TNR = TN / N
#    end
#    println("specificity, selectivity or true negative rate (TNR): $TNR")
#
#    PPV = 0
#    if (TP + FP) != 0
#        PPV = TP / (TP + FP)
#    end 
#    println("precision or positive predictive value (PPV): $PPV")
#
#    NPV = 0
#    if (TN + FN) != 0
#       NPV = TN / (TN + FN)
#    end
#    println("negative predictive value (NPV): $NPV")
#
#    FNR = 1 - TPR
#    println("miss rate or false negative rate (FNR): $FNR")
#
#    FPR = 1 - TNR
#    println("Fall-out or false positive rate (FPR): $FPR")
#
#    FDR = 1 - PPV
#    println("false discovery rate: $FDR")
#
#    FOR = 1 - NPV
#    println("false omission rate: $FOR")
#
#    ACC = (TP + TN) / total
#    println("accuracy: $ACC\n")
#
#    Fone = 0
#    if (2 * TP + FP + FN) != 0
#        Fone = (2 * TP) / (2 * TP + FP + FN)
#    end
#    println("F1 Score: $Fone")
#   
#    MCC = 0
#    divideBy = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))^(1/2)
#    if divideBy != 0
#        MCC = ((TP * TN) - (FP * FN)) / divideBy
#    end
#    println("Mathews correlation coefficient (MCC): $MCC")
#
#    BM = TPR + TNR -1
#    println("Informedness or Bookmaker Informedness (BM): $BM")
#
#    MK = PPV + NPV -1
#    println("Markedness (MK): $MK\n")
#
#    println("Expected 0 predicted 0: $expected_zero_got_zero")
#    println("Expected 0 predicted 1: $expected_zero_got_one")
#    println("Expected 0 predicted 2: $expected_zero_got_two")
#    println("Expected 1 predicted 0: $expected_one_got_zero")
#    println("Expected 1 predicted 1: $expected_one_got_one")
#    println("Expected 1 predicted 2: $expected_one_got_two")
#    println("Expected 2 predicted 0: $expected_two_got_zero")
#    println("Expected 2 predicted 1: $expected_two_got_one")
#    println("Expected 2 predicted 2: $expected_two_got_two\n\n")
#
#
#    println("Patient\tNode\tExpected\tActual\n")
#    for patientname in keys(errors)
#        patient = errors[patientname]
#        for genename in keys(patient)
#            values = errors[patientname][genename]
#            println("$patientname\t$genename\t", values[1], "\t", values[2])
#        end
#    end
end

end
