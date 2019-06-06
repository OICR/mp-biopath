module AnalyzeResultsFull


using CSV
using DataFrames
using Query
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
   push!(df, row)
end

function get_TPR(TP, P)
  if P != 0
     return TP / P
  end
  return 0
end

function get_TNR(TN, N)
  if N != 0
     return TN / N
  end
  return 0
end

function get_PPV(TP, FP)
   total = TP + FP
   if total != 0
      return TP / total
   end
   return 0
end

function get_NPV(TN, FN)
   total = TN + FN
   if total != 0
       return TN / total
   end
   return 0
end

function get_FNR(TP, P)
   return 1 - get_TPR(TP, P)
end

function get_FPR(TN, N)
   return 1 - get_TNR(TN, N)
end

function get_FDR(TP, FP)
   return 1 - get_PPV(TP, FP)
end

function get_FOR(TN, FN)
   return 1 - get_NPV(TN, FN)
end

function get_ACC(TP, TN, total)
   return (TP + TN) / total
end

function get_F1(TPR, PPV)
   if (PPV + TPR) != 0
       return 2 * ((PPV * TPR) / (PPV + TPR))
   end
   return 0
end

function addThreeDimensionalConfusionMatrixToDataFrame(df, Lowerbound, Upperbound, confusion_matrix)
   TP_down = confusion_matrix[1,1]
   FP_down = confusion_matrix[1,2] + confusion_matrix[1,3]
   FN_down = confusion_matrix[2,1] + confusion_matrix[3,1]
   TN_down = confusion_matrix[2,2] + confusion_matrix[2,3] + confusion_matrix[3,2] + confusion_matrix[3,3]

   Total = TN_down + FP_down + FN_down + TP_down
   P_down = TP_down + FN_down
   N_down = TN_down + FP_down

   TPR_down = get_TPR(TP_down, P_down)
   TNR_down = get_TNR(TN_down, N_down)
   PPV_down = get_PPV(TP_down, FP_down)
   NPV_down = get_NPV(TN_down, FN_down)
   FNR_down = get_FNR(TP_down, P_down)
   FPR_down = get_FNR(TN_down, N_down)
   FDR_down = get_FDR(TP_down, FP_down)
   FOR_down = get_FOR(TN_down, FN_down)
   ACC_down = get_ACC(TP_down, TN_down, Total)
   F1_down = get_F1(TPR_down, PPV_down)

   TP_normal = confusion_matrix[2,2]
   FP_normal = confusion_matrix[2,1] + confusion_matrix[2,3]
   FN_normal = confusion_matrix[1,2] + confusion_matrix[3,2]
   TN_normal = confusion_matrix[1,1] + confusion_matrix[1,3] + confusion_matrix[3,1] + confusion_matrix[3,3]
   P_normal = TP_normal + FN_normal
   N_normal = TN_normal + FP_normal

   TPR_normal = get_TPR(TP_normal, P_normal)
   TNR_normal = get_TNR(TN_normal, N_normal)
   PPV_normal = get_PPV(TP_normal, FP_normal)
   NPV_normal = get_NPV(TN_normal, FN_normal)
   FNR_normal = get_FNR(TP_normal, P_normal)
   FPR_normal = get_FNR(TN_normal, N_normal)
   FDR_normal = get_FDR(TP_normal, FP_normal)
   FOR_normal = get_FOR(TN_normal, FN_normal)
   ACC_normal = get_ACC(TP_normal, TN_normal, Total)
   F1_normal = get_F1(TPR_normal, PPV_normal)

   TP_up = confusion_matrix[3,3]
   FP_up = confusion_matrix[3,1] + confusion_matrix[3,2]
   FN_up = confusion_matrix[1,3] + confusion_matrix[2,3]
   TN_up = confusion_matrix[1,1] + confusion_matrix[1,2] + confusion_matrix[2,1] + confusion_matrix[2,2]
   P_up = TP_up + FN_up
   N_up = TN_up + FP_up

   TPR_up = get_TPR(TP_up, P_up)
   TNR_up = get_TNR(TN_up, N_up)
   PPV_up = get_PPV(TP_up, FP_up)
   NPV_up = get_NPV(TN_up, FN_up)
   FNR_up = get_FNR(TP_up, P_up)
   FPR_up = get_FNR(TN_up, N_up)
   FDR_up = get_FDR(TP_up, FP_up)
   FOR_up = get_FOR(TN_up, FN_up)
   ACC_up = get_ACC(TP_up, TN_up, Total)
   F1_up = get_F1(TPR_up, PPV_up)

   overall_ACC = (confusion_matrix[1,1] + confusion_matrix[2,2] + confusion_matrix[3,3]) / Total
   average_F1 = (F1_down + F1_normal + F1_up) / 3
   row = [Lowerbound, Upperbound,
          confusion_matrix[1,1], confusion_matrix[1,2], confusion_matrix[1,3],
          confusion_matrix[2,1], confusion_matrix[2,2], confusion_matrix[2,3],
          confusion_matrix[3,1], confusion_matrix[3,2], confusion_matrix[3,3],
          TP_down, FP_down, FN_down, TN_down, P_down, N_down, Total,
          TP_normal, FP_normal, FN_normal, TN_normal, P_normal, N_normal,
          TP_up, FP_up, FN_up, TN_up,
          TPR_down, TNR_down, PPV_down, NPV_down, FNR_down, FPR_down, FDR_down, FOR_down, ACC_down, F1_down,
          TPR_normal, TNR_normal, PPV_normal, NPV_normal, FNR_normal, FPR_normal, FDR_normal, FOR_normal, ACC_normal, F1_normal,
          TPR_up, TNR_up, PPV_up, NPV_up, FNR_up, FPR_up, FDR_up, FOR_up, ACC_up, F1_up,
	  overall_ACC,
          average_F1]
                    
   push!(df, row)
end

function initializeThreeDimensionalConfusionMatricies()
   return DataFrame(Lowerbound = Float64[],
                    Upperbound = Float64[],
		    ExpectedDownActualDown = Int16[],
		    ExpectedDownActualNormal = Int16[],
		    ExpectedDownActualUp = Int16[],
		    ExpectedNormalActualDown = Int16[],
		    ExpectedNormalActualNormal = Int16[],
		    ExpectedNormalActualUp = Int16[],
		    ExpectedUpActualDown = Int16[],
		    ExpectedUpActualNormal = Int16[],
		    ExpectedUpActualUp = Int16[],
                    TP_down = Int16[],
                    FP_down = Int16[],
                    FN_down = Int16[],
                    TN_down = Int16[],
                    P_down = Int16[],
                    N_down = Int16[],
                    Total = Int16[],
                    TP_normal = Int16[],
                    FP_normal = Int16[],
                    FN_normal = Int16[],
                    TN_normal = Int16[],
                    P_normal = Int16[],
                    N_normal = Int16[],
                    TP_up = Int16[],
                    FP_up = Int16[],
                    FN_up = Int16[],
                    TN_up = Int16[],
                    TPR_down = Float64[],
                    TNR_down = Float64[],
                    PPV_down = Float64[],
                    NPV_down = Float64[],
                    FNR_down = Float64[],
                    FPR_down = Float64[],
                    FDR_down = Float64[],
                    FOR_down = Float64[],
                    ACC_down = Float64[],
                    F1_down = Float64[],
                    TPR_normal = Float64[],
                    TNR_normal = Float64[],
                    PPV_normal = Float64[],
                    NPV_normal = Float64[],
                    FNR_normal = Float64[],
                    FPR_normal = Float64[],
                    FDR_normal = Float64[],
                    FOR_normal = Float64[],
                    ACC_normal = Float64[],
                    F1_normal = Float64[],
                    TPR_up = Float64[],
                    TNR_up = Float64[],
                    PPV_up = Float64[],
                    NPV_up = Float64[],
                    FNR_up = Float64[],
                    FPR_up = Float64[],
                    FDR_up = Float64[],
                    FOR_up = Float64[],
                    ACC_up = Float64[],
                    F1_up = Float64[],
		    Overall_ACC = Float64[],
                    average_F1 = Float64[])
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

function weightConfusionMatrix(confusion_matrix, total_down_expected, total_normal_expected, total_up_expected)
    total_down = confusion_matrix[1,1] + confusion_matrix[1,2] + confusion_matrix[1,3]
    total_normal = confusion_matrix[2,1] + confusion_matrix[2,2] + confusion_matrix[2,3]
    total_up = confusion_matrix[3,1] + confusion_matrix[3,2] + confusion_matrix[3,3]

    confusion_matrix_weighted = zeros(Float32, (4, 4))
    confusion_matrix_weighted[1,1] = floor(Int, confusion_matrix[1,1] / total_down * total_down_expected)
    confusion_matrix_weighted[1,2] = floor(Int, confusion_matrix[1,2] / total_down * total_down_expected)
    confusion_matrix_weighted[1,3] = floor(Int, confusion_matrix[1,3] / total_down * total_down_expected)
    confusion_matrix_weighted[2,1] = floor(Int, confusion_matrix[2,1] / total_normal * total_normal_expected)
    confusion_matrix_weighted[2,2] = floor(Int, confusion_matrix[2,2] / total_normal * total_normal_expected)
    confusion_matrix_weighted[2,3] = floor(Int, confusion_matrix[2,3] / total_normal * total_normal_expected)
    confusion_matrix_weighted[3,1] = floor(Int, confusion_matrix[3,1] / total_up * total_up_expected)
    confusion_matrix_weighted[3,2] = floor(Int, confusion_matrix[3,2] / total_up * total_up_expected)
    confusion_matrix_weighted[3,3] = floor(Int, confusion_matrix[3,3] / total_up * total_up_expected)
    return confusion_matrix_weighted
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

    total_down_expected = 0
    total_normal_expected = 0
    total_up_expected = 0
    for pathway_id in keys(expected_results)
        for scenario in keys(expected_results[pathway_id]["samplenodestate"])
            for keyoutput in keys(expected_results[pathway_id]["samplenodestate"][scenario])
                value = expected_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                if value == "0" 
                    total_down_expected += 1
                elseif value == "1"
                    total_normal_expected += 1
                elseif value == "2"
                    total_up_expected += 1
                else
                    println(value)
                    println(typeof(value))
                    exit()
                end
            end
        end
    end

    total_tests = total_down_expected + total_normal_expected + total_up_expected
    confusion_matricies = initializeThreeDimensionalConfusionMatricies()
    confusion_matricies_weighted = initializeThreeDimensionalConfusionMatricies()
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
    
    addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies, 1, 1, confusion_matrix)
    CSV.write(joinpath(analysis_results_folder, "expectedVsExperimentalConfusionMatrix.tsv"), confusion_matricies, delim = '\t')

    confusion_matrix_weighted = weightConfusionMatrix(confusion_matrix, total_down_expected, total_normal_expected, total_up_expected)
    addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies_weighted, 1, 1, confusion_matrix_weighted)
    CSV.write(joinpath(analysis_results_folder, "expectedVsExperimentalConfusionMatrixWeighted.tsv"), confusion_matricies_weighted, delim = '\t')
    #=
    ### determine optimal lowerbound cutoff based on weighted average F1 score
    confusion_matricies = initializeThreeDimensionalConfusionMatricies()
    confusion_matricies_weighted = initializeThreeDimensionalConfusionMatricies()
    for lowerbound in range(0,1, step=0.01)
        for upperbound in range(1,10, step=0.01)
            confusion_matrix = zeros(UInt16, (3, 3))
            for pathway_id in keys(experimental_results)
                for scenario in keys(experimental_results[pathway_id]["samplenodestate"])
                    for keyoutput in keys(experimental_results[pathway_id]["samplenodestate"][scenario])
                        experimental_value = experimental_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                        result_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
                        if experimental_value == "0"
                            row = 1
                        elseif experimental_value == "1"
                            row = 2
                        else
                            row = 3
                        end
                        if result_value < lowerbound
                            confusion_matrix[row,1] = confusion_matrix[row,1] + 1
                        elseif result_value > upperbound
                            confusion_matrix[row,3] = confusion_matrix[row,3] + 1
                        else
                            confusion_matrix[row,2] = confusion_matrix[row,2] + 1
                        end
                    end
                end
            end

            addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies, lowerbound, upperbound, confusion_matrix)
            confusion_matrix_weighted = weightConfusionMatrix(confusion_matrix, total_down_expected, total_normal_expected, total_up_expected)
            addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies_weighted, lowerbound, upperbound, confusion_matrix_weighted)
        end
    end

    down_confusion_matricies_weighted = @from i in confusion_matricies_weighted begin
            @where i.Upperbound == 1.05
            @select {i.Lowerbound,
                     TP=i.TP_down, FP=i.FP_down, FN=i.FN_down, TN=i.TN_down, P=i.P_down, N=i.N_down, i.Total,
                     TPR=i.TPR_down, TNR=i.TNR_down, PPV=i.PPV_down, NPV=i.NPV_down, FNR=i.FNR_down, FPR=i.FPR_down, FDR=i.FDR_down, FOR=i.FOR_down, ACC=i.ACC_down, F1=i.F1_down}
            @collect DataFrame
       end

    sort!(down_confusion_matricies_weighted, (:F1), rev=(true))
    lowerbound = down_confusion_matricies_weighted[1,:Lowerbound]

    CSV.write(joinpath(analysis_results_folder, "experimentalVsPredictedWeightedConfusionMatrix-10-pathways-down-cutoff.tsv"), down_confusion_matricies_weighted, delim = '\t')

    up_confusion_matricies_weighted = @from i in confusion_matricies_weighted begin
            @where i.Lowerbound == 0.95
            @select {i.Upperbound,
                     TP=i.TP_down, FP=i.FP_down, FN=i.FN_down, TN=i.TN_down, P=i.P_down, N=i.N_down, i.Total,
                     TPR=i.TPR_down, TNR=i.TNR_down, PPV=i.PPV_down, NPV=i.NPV_down, FNR=i.FNR_down, FPR=i.FPR_down, FDR=i.FDR_down, FOR=i.FOR_down, ACC=i.ACC_down, F1=i.F1_down}
            @collect DataFrame
       end
    sort!(up_confusion_matricies_weighted, (:F1), rev=(true))
    CSV.write(joinpath(analysis_results_folder, "experimentalVsPredictedWeightedConfusionMatrix-10-pathways-down-cutoff.tsv"), up_confusion_matricies_weighted, delim = '\t')
    upperbound = up_confusion_matricies_weighted[1,:Upperbound]

    CSV.write(joinpath(analysis_results_folder, "experimentalVsPredictedWeightedConfusionMatrix-10-pathways.tsv"), confusion_matricies_weighted, delim = '\t')
    CSV.write(joinpath(analysis_results_folder, "experimentalVsPredictedConfusionMatrix-10-pathways.tsv"), confusion_matricies, delim = '\t')

    println("Ideal lowerbound calculated to be: $lowerbound")
    println("Ideal upperbound calculated to be: $upperbound")

    optimal_confusion_matrix_weighted = @from i in confusion_matricies_weighted begin
            @where i.Lowerbound == lowerbound && i.Upperbound == upperbound
            @select i
            @collect DataFrame
    end

    #create based on experimental pathways
    confusion_matricies = initializeThreeDimensionalConfusionMatricies()
    confusion_matrix = zeros(UInt16, (3, 3))
    for pathway_id in keys(experimental_results)
        for scenario in keys(experimental_results[pathway_id]["samplenodestate"])
            for keyoutput in keys(experimental_results[pathway_id]["samplenodestate"][scenario])
                expected_value = expected_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                result_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
                if expected_value == "0"
                    row = 1
                elseif expected_value == "1"
                    row = 2
                else
                    row = 3
                end
                if result_value < lowerbound
                    confusion_matrix[row,1] = confusion_matrix[row,1] + 1
                elseif result_value > upperbound
                    confusion_matrix[row,3] = confusion_matrix[row,3] + 1
                else
                    confusion_matrix[row,2] = confusion_matrix[row,2] + 1
                end
            end
        end
    end
    addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies, lowerbound, upperbound, confusion_matrix)
    CSV.write(joinpath(analysis_results_folder, "expectedVsPredictedConfusionMatrix-10-pathways.tsv"), confusion_matricies, delim = '\t')
    =#
    lowerbound = 0.83
    upperbound = 1.08
    # create based on all patwhays
    confusion_matricies = initializeThreeDimensionalConfusionMatricies()
    confusion_matrix = zeros(UInt16, (3, 3))
    for pathway_id in keys(expected_results)
        for scenario in keys(expected_results[pathway_id]["samplenodestate"])
            for keyoutput in keys(expected_results[pathway_id]["samplenodestate"][scenario])
                expected_value = expected_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                result_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
                if expected_value == "0"
                    row = 1
                elseif expected_value == "1"
                    row = 2
                else
                    row = 3
                end
                if result_value < lowerbound
                    confusion_matrix[row,1] = confusion_matrix[row,1] + 1
                elseif result_value > upperbound
                    confusion_matrix[row,3] = confusion_matrix[row,3] + 1
                else
                    confusion_matrix[row,2] = confusion_matrix[row,2] + 1
                end
            end
        end
    end
    addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies, lowerbound, upperbound, confusion_matrix)
    CSV.write(joinpath(analysis_results_folder, "expectedVsPredictedConfusionMatrix-all-pathways.tsv"), confusion_matricies, delim = '\t')

    pathway_tests = CSV.read("/data/MP_BioPathReactomePathwayAccuracy-Tests.tsv", delim="\t",copycols=true)
    for pathway_id in keys(experimental_results)
        for scenario in keys(expected_results[pathway_id]["samplenodestate"])
            for keyoutput in keys(expected_results[pathway_id]["samplenodestate"][scenario])
                #expected_value = expected_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                result_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
                result_state = "normal"
                if result_value < lowerbound
                    result_state = "down"
                elseif result_value > upperbound
                    result_state = "up"
                end
                pathway_tests[((pathway_tests[:scenario] .== scenario) .& (pathway_tests[:keyoutput_id] .== parse(Int64, keyoutput)) .& (pathway_tests[:pathway_id] .== parse(Int64, pathway_id))),:mp_biopath_state] = result_state
                pathway_tests[((pathway_tests[:scenario] .== scenario) .& (pathway_tests[:keyoutput_id] .== parse(Int64, keyoutput)) .& (pathway_tests[:pathway_id] .== parse(Int64, pathway_id))),:mp_biopath_value] = result_value
            end
        end
    end
    CSV.write(joinpath(analysis_results_folder, "MP_BioPathReactomePathwayAccuracy-Tests-updated.tsv"), pathway_tests, delim = '\t')
end

end
