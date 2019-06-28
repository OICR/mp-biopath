module AnalyzeResultsFull


using CSV
using DataFrames
using Query
using Statistics
include("valuetostate.jl")
include("probability.jl")
include("evidence.jl")
include("results.jl")
include("idMap.jl")
include("pathwayList.jl")
include("pi.jl")

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
   weighted_F1 = (F1_down * P_down + F1_normal * P_normal + F1_up * P_up) / Total
   row = [Lowerbound, Upperbound,
          confusion_matrix[1,1], confusion_matrix[1,2], confusion_matrix[1,3],
          confusion_matrix[2,1], confusion_matrix[2,2], confusion_matrix[2,3],
          confusion_matrix[3,1], confusion_matrix[3,2], confusion_matrix[3,3],
          TP_down, FP_down, FN_down, TN_down, P_down, N_down, Total,
          TP_normal, FP_normal, FN_normal, TN_normal, P_normal, N_normal,
          TP_up, FP_up, FN_up, TN_up, P_up, N_up,
          TPR_down, TNR_down, PPV_down, NPV_down, FNR_down, FPR_down, FDR_down, FOR_down, ACC_down, F1_down,
          TPR_normal, TNR_normal, PPV_normal, NPV_normal, FNR_normal, FPR_normal, FDR_normal, FOR_normal, ACC_normal, F1_normal,
          TPR_up, TNR_up, PPV_up, NPV_up, FNR_up, FPR_up, FDR_up, FOR_up, ACC_up, F1_up,
	  overall_ACC,
          average_F1, weighted_F1]
                    
   push!(df, row)
end

function initializeThreeDimensionalConfusionMatricies()
   return DataFrame(Lowerbound = Float64[],
                    Upperbound = Float64[],
		    PredictedDownActualDown = Int16[],
		    PredictedDownActualNormal = Int16[],
		    PredictedDownActualUp = Int16[],
		    PredictedNormalActualDown = Int16[],
		    PredictedNormalActualNormal = Int16[],
		    PredictedNormalActualUp = Int16[],
		    PredictedUpActualDown = Int16[],
		    PredictedUpActualNormal = Int16[],
		    PredictedUpActualUp = Int16[],
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
		    P_up = Int16[],
		    N_up = Int16[],
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
                    average_F1 = Float64[],
                    weighted_F1 = Float64[])
end

function weightConfusionMatrix(confusion_matrix, total_down_reactome_curator, total_normal_reactome_curator, total_up_reactome_curator)
    total_down = confusion_matrix[1,1] + confusion_matrix[1,2] + confusion_matrix[1,3]
    total_normal = confusion_matrix[2,1] + confusion_matrix[2,2] + confusion_matrix[2,3]
    total_up = confusion_matrix[3,1] + confusion_matrix[3,2] + confusion_matrix[3,3]

    confusion_matrix_weighted = zeros(Float32, (4, 4))
    confusion_matrix_weighted[1,1] = floor(Int, confusion_matrix[1,1] / total_down * total_down_reactome_curator)
    confusion_matrix_weighted[1,2] = floor(Int, confusion_matrix[1,2] / total_down * total_down_reactome_curator)
    confusion_matrix_weighted[1,3] = floor(Int, confusion_matrix[1,3] / total_down * total_down_reactome_curator)
    confusion_matrix_weighted[2,1] = floor(Int, confusion_matrix[2,1] / total_normal * total_normal_reactome_curator)
    confusion_matrix_weighted[2,2] = floor(Int, confusion_matrix[2,2] / total_normal * total_normal_reactome_curator)
    confusion_matrix_weighted[2,3] = floor(Int, confusion_matrix[2,3] / total_normal * total_normal_reactome_curator)
    confusion_matrix_weighted[3,1] = floor(Int, confusion_matrix[3,1] / total_up * total_up_reactome_curator)
    confusion_matrix_weighted[3,2] = floor(Int, confusion_matrix[3,2] / total_up * total_up_reactome_curator)
    confusion_matrix_weighted[3,3] = floor(Int, confusion_matrix[3,3] / total_up * total_up_reactome_curator)
    return confusion_matrix_weighted
end

function createConfusionMatrixFormatted(cm, actual_class, predicted_class, filepath)
    predicted_down_total = cm[:PredictedDownActualDown] + cm[:PredictedDownActualNormal] + cm[:PredictedDownActualUp] 
    predicted_normal_total = cm[:PredictedNormalActualDown] + cm[:PredictedNormalActualNormal] + cm[:PredictedNormalActualUp] 
    predicted_up_total = cm[:PredictedUpActualDown] + cm[:PredictedUpActualNormal] + cm[:PredictedUpActualUp] 
    actual_down_total = cm[:PredictedDownActualDown] + cm[:PredictedNormalActualDown] + cm[:PredictedUpActualDown] 
    actual_normal_total = cm[:PredictedDownActualNormal] + cm[:PredictedNormalActualNormal] + cm[:PredictedUpActualNormal] 
    actual_up_total = cm[:PredictedDownActualUp] + cm[:PredictedNormalActualUp] + cm[:PredictedUpActualUp] 
    title =  string(predicted_class, " vs ", actual_class)
    A = [title "" "" actual_class "" "" "Precision" "";
         "" "" "Down" "Normal" "Up" "Total" "(PPV)" "F1";
         "" "Down" cm[:PredictedDownActualDown] cm[:PredictedDownActualNormal] cm[:PredictedDownActualUp]  predicted_down_total cm[:PPV_down] cm[:F1_down];
         predicted_class "Normal" cm[:PredictedNormalActualDown] cm[:PredictedNormalActualNormal] cm[:PredictedNormalActualUp]  predicted_normal_total cm[:PPV_normal] cm[:F1_normal];
         "" "Up" cm[:PredictedUpActualDown] cm[:PredictedUpActualNormal] cm[:PredictedUpActualUp] predicted_up_total cm[:PPV_up] cm[:F1_up];
         "" "Total" actual_down_total actual_normal_total actual_up_total cm[:Total] "" "";
         "Sensitivity" "(TPR)" cm[:TPR_down] cm[:TPR_normal] cm[:TPR_up] "" cm[:Overall_ACC] cm[:average_F1];
         "" "" "" "" "" "" "" "";
         "Predicted Down" "" "Positive" "Negative" "ACC" "" "" "";
         "" "Positive" cm[:TP_down] cm[:FP_down] cm[:ACC_down] "" "" "";
         "" "Negative" cm[:FN_down] cm[:TN_down] "" "" "" "";
         "" "" "" "" "" "" "" "";
         "Predicted Normal" "" "Positive" "Negative" "ACC" "" "" "";
         "" "Positive" cm[:TP_normal] cm[:FP_normal] cm[:ACC_normal] "" "" "";
         "" "Negative" cm[:FN_normal] cm[:TN_normal] "" "" "" "";
         "" "" "" "" "" "" "" "";
         "Predicted Up" "" "Positive" "Negative" "ACC" "" "" "";
         "" "Positive" cm[:TP_up] cm[:FP_up] cm[:ACC_up] "" "" "";
         "" "Negative" cm[:FN_up] cm[:TN_up] "" "" "" ""
         "" "" "" "" "" "" "" "";
         "" "TP" "TN" "FP" "FN" "Total" "ACC" "";
         "Up"  cm[:TP_up] cm[:TN_up] cm[:FP_up] cm[:FN_up] cm[:Total] cm[:ACC_up] "";
         "Normal"  cm[:TP_normal] cm[:TN_normal] cm[:FP_normal] cm[:FN_normal] cm[:Total] cm[:ACC_normal] "";
         "Down" cm[:TP_down] cm[:TN_down] cm[:FP_down] cm[:FN_down] cm[:Total] cm[:ACC_down] ""]

    CSV.write(filepath,  DataFrame(A), writeheader=false)
 end

function getParents(pathway_network, node_id)
    node = pathway_network[string(node_id)]

    parents = []
    if node.relation in ["AND", "NEG", "ANDNEG"]
        parents = node.parents
        node.parents
    elseif node.relation == "OR"
        parents = node.posParents
    end
    
    non_pseudo_parents = Set()
    for parent_id in parents
       if occursin("PSEUDONODE", parent_id) == true
           non_pseudo_parents_of_parents = getParents(pathway_network, parent_id)
           union!(non_pseudo_parents, non_pseudo_parents_of_parents)
       else
           push!(non_pseudo_parents, parent_id)
       end
    end
    return non_pseudo_parents
end

function findDistanceToInput(pathway_network, output_node_id, input_id)
    distance = 0
    first_iteration = true
    parents = Set()
    visited = Set([input_id])
    while first_iteration == true || length(parents) != 0
        if output_node_id == input_id
            return distance
        end
        if first_iteration == true
            parents = Set(getParents(pathway_network, output_node_id))
        end
        distance += 1
        new_parents = Set()
        for parent_id in parents
             if parent_id == input_id
                 return distance
             elseif in(parent_id, visited) == false
                 push!(visited, parent_id)
                 union!(new_parents, getParents(pathway_network, parent_id))
             end
        end
               
        parents = new_parents
        first_iteration = false
    end
    return -1
end

function analyzeResultsFull(results_folder, reactome_curator_results_folder, pathway_list_file, pathways_folder, db_id_name_mapping_file, analysis_results_folder, experimental_results_folder, verbose)
    if isfile(pathway_list_file)
        pathway_name_to_id_map = PathwayList.getPathwayNameToIdMap(pathway_list_file)
    else
        println("ERROR pathway list file does not exists: $pathway_list")
        exit(2)
    end

    if isfile(db_id_name_mapping_file)
        id_map = IdMap.getIDmap(db_id_name_mapping_file)
    else
        println("ERROR db_id_name_mapping_file does not exist: $db_id_name_mapping_file")
        exit(3)
    end
    reactome_curator_results_file_ending = "_reactome_curator_results.tsv"
    reactome_curator_results = Dict()
    if isdir(reactome_curator_results_folder)
        for reactome_curator_results_filename in readdir(reactome_curator_results_folder)
            if endswith(reactome_curator_results_filename, reactome_curator_results_file_ending)
                reactome_curator_results_filepath = joinpath(reactome_curator_results_folder, reactome_curator_results_filename)
                pathway_reactome_curator_results = Results.getExpected(reactome_curator_results_filepath)
                pathway_name = reactome_curator_results_filename[1:end - lastindex(reactome_curator_results_file_ending)]
                pathway_id = pathway_name_to_id_map[pathway_name]
                reactome_curator_results[pathway_id] = pathway_reactome_curator_results
             end
        end
    else
        println("ERROR: Folder does not exist $reactome_curator_results_folder")
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
 
    pathway_tests = CSV.read("/data/MP_BioPathReactomePathwayAccuracy-Tests.tsv", delim="\t",copycols=true)
    input_to_output_distance = DataFrame(pathway_id = String[],
	                                 input_name = String[],
					 input_id = String[],
					 keyoutput_id = String[],
					 distance = String[])
    for i in 1:size(pathway_tests,1)
        pathway_id = pathway_tests[i,:pathway_id]
        pathway_network = pathway_networks[string(pathway_id)]
        scenario = pathway_tests[i,:scenario]
        keyoutput_id = pathway_tests[i,:keyoutput_id]
        experimental_value = pathway_tests[i,:experimental_value]

        expected_value = pathway_tests[i,:expected_value]

        reactome_curator_result = reactome_curator_results[string(pathway_id)]["samplenodestate"][scenario][string(keyoutput_id)]

        if reactome_curator_result != string(expected_value)
            println("Curator expected value does not match for pathway $pathway_id, scenario $scenario, keyoutput_id $keyoutput_id: reactome_curator: $reactome_curator_result, from_old_file: $expected_value")
        end

        if experimental_value != "NA"
            experimental_result = experimental_results[string(pathway_id)]["samplenodestate"][scenario][string(keyoutput_id)]
            if experimental_result != experimental_value
                println("Experimental value does not match for pathway $pathway_id, scenario $scenario, keyoutput_id $keyoutput_id: experimental value: $experimental_result, from_old_file: $experimental_value")
            end
        end
        input_name = scenario[1:end-2]
        for id_map_row in id_map[input_name]
            input_id = id_map_row[:Database_Identifier]
            if input_id in keys(pathway_network) && pathway_network[string(input_id)].relation == "ROOT"
                distance_rows = @from i in input_to_output_distance begin
                                    @where i.pathway_id == string(pathway_id) && i.input_name == string(input_name) && i.keyoutput_id == string(keyoutput_id)
                                    @select i
                                    @collect DataFrame
                                end
                if (length(distance_rows[:,1]) == 0)
                    distance = findDistanceToInput(pathway_network, keyoutput_id, input_id)
                    if distance != -1
                        distance_str = string(distance)
                        push!(input_to_output_distance, [string(pathway_id), string(input_name), string(input_id), string(keyoutput_id), distance_str])
                    end
		end
	    end
        end
    end
    CSV.write(joinpath(analysis_results_folder, "input_to_keyoutput_distance.tsv"), input_to_output_distance, delim = '\t')

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

    total_down_reactome_curator = 0
    total_normal_reactome_curator = 0
    total_up_reactome_curator = 0
    for pathway_id in keys(reactome_curator_results)
        for scenario in keys(reactome_curator_results[pathway_id]["samplenodestate"])
            if scenario == "control"
                continue
            end
            for keyoutput in keys(reactome_curator_results[pathway_id]["samplenodestate"][scenario])
                value = reactome_curator_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                if value == "0" 
                    total_down_reactome_curator += 1
                elseif value == "1"
                    total_normal_reactome_curator += 1
                elseif value == "2"
                    total_up_reactome_curator += 1
                else
                    println(value)
                    println(typeof(value))
                    exit()
                end
            end
        end
    end

    total_tests = total_down_reactome_curator + total_normal_reactome_curator + total_up_reactome_curator
#    confusion_matricies = initializeThreeDimensionalConfusionMatricies()
#    confusion_matricies_weighted = initializeThreeDimensionalConfusionMatricies()
#    confusion_matrix = zeros(Int16, (4, 4))
#    for pathway_id in keys(experimental_results)
#        for scenario in keys(experimental_results[pathway_id]["samplenodestate"])
#            for keyoutput in keys(experimental_results[pathway_id]["samplenodestate"][scenario])
#                  reactome_curator_value = reactome_curator_results[pathway_id]["samplenodestate"][scenario][keyoutput]
#                  experimental_value = experimental_results[pathway_id]["samplenodestate"][scenario][keyoutput]
#                  if reactome_curator_value == "0" &&experimental_value == "0"
#                      confusion_matrix[1,1] = confusion_matrix[1,1] + 1
#                  elseif reactome_curator_value == "1" && experimental_value == "0"
#                      confusion_matrix[2,1] = confusion_matrix[2,1] + 1
#                  elseif reactome_curator_value == "2" && experimental_value == "0"
#                      confusion_matrix[3,1] = confusion_matrix[3,1] + 1
#                  elseif reactome_curator_value == "0" && experimental_value == "1"
#                      confusion_matrix[1,2] = confusion_matrix[1,2] + 1
#                  elseif reactome_curator_value == "1" && experimental_value == "1"
#                      confusion_matrix[2,2] = confusion_matrix[2,2] + 1
#                  elseif reactome_curator_value == "2" && experimental_value == "1"
#                      confusion_matrix[3,2] = confusion_matrix[3,2] + 1
#                  elseif reactome_curator_value == "0" && experimental_value == "2"
#                      confusion_matrix[1,3] = confusion_matrix[1,3] + 1
#                  elseif reactome_curator_value == "1" && experimental_value == "2"
#                      confusion_matrix[2,3] = confusion_matrix[2,3] + 1
#                  elseif reactome_curator_value == "2" && experimental_value == "2"
#                      confusion_matrix[3,3] = confusion_matrix[3,3] + 1
#                  end
#             end
#         end
#    end
#    
#    addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies, 1, 1, confusion_matrix)
#    createConfusionMatrixFormatted(confusion_matricies, "Experimental", "Reactome Curator", joinpath(analysis_results_folder, "reactomeCuratorVsExperimentalConfusionMatrix.csv"))
#
#    confusion_matrix_weighted = weightConfusionMatrix(confusion_matrix, total_down_reactome_curator, total_normal_reactome_curator, total_up_reactome_curator)
#    addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies_weighted, 1, 1, confusion_matrix_weighted)
#    createConfusionMatrixFormatted(confusion_matricies_weighted, "Experimental", "Reactome Curator", joinpath(analysis_results_folder, "reactomeCuratorVsExperimentalConfusionMatrixWeighted.csv"))
#
#    ### determine optimal lowerbound cutoff based on weighted average F1 score
#    confusion_matricies = initializeThreeDimensionalConfusionMatricies()
#    confusion_matricies_weighted = initializeThreeDimensionalConfusionMatricies()
#    for lowerbound in range(0.1,1, step=0.01)
#        for upperbound in range(1.01,10, step=0.01)
#            confusion_matrix = zeros(UInt16, (3, 3))
#            for pathway_id in keys(experimental_results)
#                for scenario in keys(experimental_results[pathway_id]["samplenodestate"])
#                    for keyoutput in keys(experimental_results[pathway_id]["samplenodestate"][scenario])
#                        experimental_value = experimental_results[pathway_id]["samplenodestate"][scenario][keyoutput]
#                        result_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
#                        if experimental_value == "0"
#                            column = 1
#                        elseif experimental_value == "1"
#                            column = 2
#                        else
#                            column = 3
#                        end
#                        if result_value < lowerbound
#                            confusion_matrix[1,column] = confusion_matrix[1,column] + 1
#                        elseif result_value > upperbound
#                            confusion_matrix[3,column] = confusion_matrix[3,column] + 1
#                        else
#                            confusion_matrix[2,column] = confusion_matrix[2,column] + 1
#                        end
#                    end
#                end
#            end
#
#            addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies, lowerbound, upperbound, confusion_matrix)
#            confusion_matrix_weighted = weightConfusionMatrix(confusion_matrix, total_down_reactome_curator, total_normal_reactome_curator, total_up_reactome_curator)
#            addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies_weighted, lowerbound, upperbound, confusion_matrix_weighted)
#        end
#    end
#
#    CSV.write(joinpath(analysis_results_folder, "experimentalVsMP-BioPathWeightedConfusionMatricies-10-pathways.tsv"), confusion_matricies_weighted, delim = '\t')
#    CSV.write(joinpath(analysis_results_folder, "experimentalVsMP-BioPathConfusionMatricies-10-pathways.tsv"), confusion_matricies, delim = '\t')
#
#    down_confusion_matricies_weighted = @from i in confusion_matricies_weighted begin
#            @where i.Upperbound == 1.05
#            @select {i.Lowerbound,
#                     TP=i.TP_down, FP=i.FP_down, FN=i.FN_down, TN=i.TN_down, P=i.P_down, N=i.N_down, i.Total,
#                     TPR=i.TPR_down, TNR=i.TNR_down, PPV=i.PPV_down, NPV=i.NPV_down, FNR=i.FNR_down, FPR=i.FPR_down, FDR=i.FDR_down, FOR=i.FOR_down, ACC=i.ACC_down, F1=i.F1_down}
#            @collect DataFrame
#       end
#
#    sort!(down_confusion_matricies_weighted, (:F1), rev=(true))
#    lowerbound = down_confusion_matricies_weighted[1,:Lowerbound]
#    CSV.write(joinpath(analysis_results_folder, "MP-BioPathWeightedVsExperimentalConfusionMatricies-10-pathways-down-cutoff.tsv"), down_confusion_matricies_weighted, delim = '\t')
#
#    up_confusion_matricies_weighted = @from i in confusion_matricies_weighted begin
#            @where i.Lowerbound == 0.95
#            @select {i.Upperbound,
#                     TP=i.TP_up, FP=i.FP_up, FN=i.FN_up, TN=i.TN_up, P=i.P_up, N=i.N_up, i.Total,
#                     TPR=i.TPR_up, TNR=i.TNR_up, PPV=i.PPV_up, NPV=i.NPV_up, FNR=i.FNR_up, FPR=i.FPR_up, FDR=i.FDR_up, FOR=i.FOR_up, ACC=i.ACC_up, F1=i.F1_up}
#            @collect DataFrame
#       end
#    sort!(up_confusion_matricies_weighted, (:F1), rev=(true))
#    CSV.write(joinpath(analysis_results_folder, "experimentalVsMP-BioPathWeightedConfusionMatricies-10-pathways-up-cutoff.tsv"), up_confusion_matricies_weighted, delim = '\t')
#    upperbound = up_confusion_matricies_weighted[1,:Upperbound]

    lowerbound = 0.92
    upperbound = 1.03

    println("Ideal lowerbound calculated to be: $lowerbound")
    println("Ideal upperbound calculated to be: $upperbound")

#    optimal_confusion_matrix = @from i in confusion_matricies begin
#            @where i.Lowerbound == lowerbound && i.Upperbound == upperbound
#            @select i
#            @collect DataFrame
#    end
#
#    createConfusionMatrixFormatted(optimal_confusion_matrix,
#				   "Experimental",
#				   "MP-BioPath",
#				   joinpath(analysis_results_folder, "MP-BioPathVsExperimentalConfusionMatrix.csv"))
#
#    optimal_confusion_matrix_weighted = @from i in confusion_matricies_weighted begin
#            @where i.Lowerbound == lowerbound && i.Upperbound == upperbound
#            @select i
#            @collect DataFrame
#    end
#
#    createConfusionMatrixFormatted(optimal_confusion_matrix_weighted,
#				   "Experimental",
#				   "MP-BioPath Weighted",
#				   joinpath(analysis_results_folder, "MP-BioPathWeightedVsExperimentalConfusionMatrix.csv"))
#
    #create based on experimental pathways
    confusion_matricies = initializeThreeDimensionalConfusionMatricies()
    confusion_matrix = zeros(UInt16, (3, 3))
    for pathway_id in keys(experimental_results)
        for scenario in keys(experimental_results[pathway_id]["samplenodestate"])
            for keyoutput in keys(experimental_results[pathway_id]["samplenodestate"][scenario])
                reactome_curator_value = reactome_curator_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                result_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
                if reactome_curator_value == "0"
                    column = 1
                elseif reactome_curator_value == "1"
                    column = 2
                else
                    column = 3
                end
                if result_value < lowerbound
                    confusion_matrix[1,column] = confusion_matrix[1,column] + 1
                elseif result_value > upperbound
                    confusion_matrix[3,column] = confusion_matrix[3,column] + 1
                else
                    confusion_matrix[2,column] = confusion_matrix[2,column] + 1
                end
            end
        end
    end
    addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies, lowerbound, upperbound, confusion_matrix)
    CSV.write(joinpath(analysis_results_folder, "reactomeCuratorVsMP-BioPathConfusionMatrix-10-pathways-with-tests-with-experimental-evidence.tsv"), confusion_matricies, delim = '\t')
    createConfusionMatrixFormatted(confusion_matricies,
				   "MP-BioPath",
				   "Reactome Curator",
				   joinpath(analysis_results_folder, "MP-BioPathVsReactomeCurator-10-pathways-with-tests-with-experimental-evidenceConfusionMatrix.csv"))

    confusion_matricies = initializeThreeDimensionalConfusionMatricies()
    confusion_matrix = zeros(UInt16, (3, 3))
    for pathway_id in keys(experimental_results)
        for scenario in keys(reactome_curator_results[pathway_id]["samplenodestate"])
            if scenario == "control"
                continue
            end
            for keyoutput in keys(reactome_curator_results[pathway_id]["samplenodestate"][scenario])
                reactome_curator_value = reactome_curator_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                result_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
                if reactome_curator_value == "0"
                    column = 1
                elseif reactome_curator_value == "1"
                    column = 2
                else
                    column = 3
                end
                if result_value < lowerbound
                    confusion_matrix[1,column] = confusion_matrix[1,column] + 1
                elseif result_value > upperbound
                    confusion_matrix[3,column] = confusion_matrix[3,column] + 1
                else
                    confusion_matrix[2,column] = confusion_matrix[2,column] + 1
                end
            end
        end
    end
    addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies, lowerbound, upperbound, confusion_matrix)
    createConfusionMatrixFormatted(confusion_matricies,
				   "Reactome Curator",
				   "MP-BioPath",
				   joinpath(analysis_results_folder, "MP-BioPathVsReactomeCurator-10-pathways-with-all-tests-ConfusionMatrix.csv"))

    # create based on all patwhays
    confusion_matricies = initializeThreeDimensionalConfusionMatricies()
    confusion_matrix = zeros(UInt16, (3, 3))
    for pathway_id in keys(reactome_curator_results)
        for scenario in keys(reactome_curator_results[pathway_id]["samplenodestate"])
            if scenario == "control"
                continue
            end
            for keyoutput in keys(reactome_curator_results[pathway_id]["samplenodestate"][scenario])
                reactome_curator_value = reactome_curator_results[pathway_id]["samplenodestate"][scenario][keyoutput]
                result_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
                if reactome_curator_value == "0"
                    column = 1
                elseif reactome_curator_value == "1"
                    column = 2
                else
                    column = 3
                end
                if result_value < lowerbound
                    confusion_matrix[1,column] = confusion_matrix[1,column] + 1
                elseif result_value > upperbound
                    confusion_matrix[3,column] = confusion_matrix[3,column] + 1
                else
                    confusion_matrix[2,column] = confusion_matrix[2,column] + 1
                end
            end
        end
    end
    addThreeDimensionalConfusionMatrixToDataFrame(confusion_matricies, lowerbound, upperbound, confusion_matrix)
    createConfusionMatrixFormatted(confusion_matricies,
				   "Reactome Curator",
				   "MP-BioPath",
				   joinpath(analysis_results_folder, "MP-BioPathVsReactomeCurator-all-pathways-ConfusionMatrix.csv"))

    pathway_tests[:computed_distance] = -1
    for pathway_id in keys(experimental_results)
        for scenario in keys(reactome_curator_results[pathway_id]["samplenodestate"])
            for keyoutput in keys(reactome_curator_results[pathway_id]["samplenodestate"][scenario])
                result_value = parse(Float64, mp_biopath_results[pathway_id][scenario][keyoutput])
                result_state = 1
                if result_value < lowerbound
                    result_state = 0
                elseif result_value > upperbound
                    result_state = 2
                end
                input_name = scenario[1:end-2]
                pathway_tests[((pathway_tests[:scenario] .== scenario) .& (pathway_tests[:keyoutput_id] .== parse(Int64, keyoutput)) .& (pathway_tests[:pathway_id] .== parse(Int64, pathway_id))),:mp_biopath_state] = result_state
                pathway_tests[((pathway_tests[:scenario] .== scenario) .& (pathway_tests[:keyoutput_id] .== parse(Int64, keyoutput)) .& (pathway_tests[:pathway_id] .== parse(Int64, pathway_id))),:mp_biopath_value] = result_value
                distance_rows = @from i in input_to_output_distance begin
                                    @where i.pathway_id == pathway_id && i.input_name == input_name && i.keyoutput_id == keyoutput
                                    @select i
                                    @collect DataFrame
                                end
                distances = []
                for x = 1:length(distance_rows[:,1])
                    distance = distance_rows[x,:distance]
                    if distance != "NA"
                        push!(distances, parse(Int64, distance))
                    end
                end
                if length(distances) == 0
                   distance = -1
                else
                   distance = mean(distances)
                end
                pathway_tests[((pathway_tests[:scenario] .== scenario) .& (pathway_tests[:keyoutput_id] .== parse(Int64, keyoutput)) .& (pathway_tests[:pathway_id] .== parse(Int64, pathway_id))),:computed_distance] = distance
            end
        end
    end
    CSV.write(joinpath(analysis_results_folder, "MP_BioPathReactomePathwayAccuracy-Tests-updated.tsv"), pathway_tests, delim = '\t')

    #Create tables for pathway_distance 
    ## MP-Biopath vs Experimental
    ## Reactome Curator vs Experimental
    expected_vs_actual_distances = DataFrame(expected_state = Int64[],
                                             actual_state = String[],
                                             count_with_path_mp_biopath = Int64[],
                                             count_no_path_mp_biopath = Int64[],
                                             mean_distance_mp_biopath = Float64[],
                                             variance_distance_mp_biopath = Float64[],
                                             count_with_path_reactome_curator = Int64[],
                                             count_no_path_reactome_curator = Int64[],
                                             mean_distance_reactome_curator = Float64[],
                                             variance_distance_reactome_curator = Float64[])
    for actual_state in ["0","1", "2", "NA"]
        for expected_state in [0, 1, 2] #&& i.experimental_value == actual && i.distance == "NA"
            reactome_curator_no_path = @from i in pathway_tests begin
                                           @where i.expected_value == expected_state && i.experimental_value == actual_state && i.computed_distance == -1
                                           @select i
                                           @collect DataFrame
                                       end
            count_no_path_reactome_curator = length(reactome_curator_no_path[:,1])
            reactome_curator_with_path = @from i in pathway_tests begin
                                             @where i.expected_value == expected_state && i.experimental_value == actual_state && i.computed_distance != -1
                                             @select i
                                             @collect DataFrame
                                         end
            count_with_path_reactome_curator = length(reactome_curator_with_path[:,1])
            mean_distance_reactome_curator = mean(reactome_curator_with_path[:computed_distance])
            variance_distance_reactome_curator = var(reactome_curator_with_path[:computed_distance])

            mp_biopath_no_path = @from i in pathway_tests begin
                                           @where i.mp_biopath_state == expected_state && i.experimental_value == actual_state && i.computed_distance == -1
                                           @select i
                                           @collect DataFrame
                                       end
            count_no_path_mp_biopath = length(mp_biopath_no_path[:,1])
            mp_biopath_with_path = @from i in pathway_tests begin
                                             @where i.mp_biopath_state == expected_state && i.experimental_value == actual_state && i.computed_distance != -1
                                             @select i
                                             @collect DataFrame
                                         end
            count_with_path_mp_biopath = length(mp_biopath_with_path[:,1])
            mean_distance_mp_biopath = mean(mp_biopath_with_path[:computed_distance])
            variance_distance_mp_biopath = var(mp_biopath_with_path[:computed_distance])
            push!(expected_vs_actual_distances, [expected_state,
			 		         actual_state,
						 count_with_path_mp_biopath,
						 count_no_path_mp_biopath,
						 mean_distance_mp_biopath,
						 variance_distance_mp_biopath,
						 count_with_path_reactome_curator,
						 count_no_path_reactome_curator,
						 mean_distance_reactome_curator,
						 variance_distance_reactome_curator])
        end
    end
    CSV.write(joinpath(analysis_results_folder, "expected_vs_actual_distances.tsv"), expected_vs_actual_distances, delim = '\t')
end

end
