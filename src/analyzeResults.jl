module AnalyzeResults

include("valuetostate.jl")
include("probability.jl")
include("evidence.jl")
include("results.jl")
include("idMap.jl")

function analyzeResults(resultsfile, expectedfile, downregulatedcutoff, upregulatedcutoff, verbose)
    expected_data = Results.getExpected(expectedfile)
    expected = expected_data["samplenodestate"]

    expected_zero_got_zero = 0
    expected_one_got_one = 0
    expected_two_got_two = 0
    expected_zero_got_one = 0
    expected_zero_got_two = 0
    expected_one_got_zero = 0
    expected_one_got_two = 0
    expected_two_got_one = 0
    expected_two_got_zero = 0

    results = Results.getResults(resultsfile)
    probability = 1;
    errors = Dict()
    for patientname in keys(expected)
        expected_patient_nodes = expected[patientname]
        results_patient_nodes = results[patientname]
        for nodename in keys(expected_patient_nodes)
            expectedState = expected_patient_nodes[nodename]
            resultValue = results_patient_nodes[nodename]

            resultState = ValueToState.getStateNumber(resultValue, downregulatedcutoff, upregulatedcutoff)
            if expectedState != resultState
                if haskey(errors, patientname) == false
                    errors[patientname] = Dict()
                end
                errors[patientname][nodename] = [expectedState, resultValue]
            end

            if expectedState == "0" && resultState == "0"
                expected_zero_got_zero += 1
            elseif expectedState == "1" && resultState == "1"
                expected_one_got_one += 1
            elseif expectedState == "2" && resultState == "2"
                expected_two_got_two += 1
            elseif expectedState == "0" && resultState == "1"
                expected_zero_got_one += 1
            elseif expectedState == "0" && resultState == "2"
                expected_zero_got_two += 1
            elseif expectedState == "1" && resultState == "0"
                expected_one_got_zero += 1
            elseif expectedState == "1" && resultState == "2"
                expected_one_got_two += 1
            elseif expectedState == "2" && resultState == "1"
                expected_two_got_one += 1
            elseif expectedState == "2" && resultState == "0"
                expected_two_got_zero += 1
            end
        end
    end

    TN = expected_one_got_one
    println("\ntrue negative (TN): $TN")

    TP = expected_two_got_two + expected_zero_got_zero
    println("true positive (TP): $TP")

    FP = expected_one_got_zero + expected_one_got_two + expected_two_got_zero + expected_zero_got_two
    println("false positive (FP): $FP")

    FN = expected_two_got_one + expected_zero_got_one
    println("false negative (FN): $FN")

    total = TN + FP + FN + TP
    println("Total number of cases: $total")

    P = TP + FN 
    println("conditional positive (P): $P")    

    N = TN + FP
    println("conditional negative (N): $N")

    TPR = 0
    if P != 0
        TPR = TP / P
    end
    println("sensitivity, recall, hit rate, or true positive rate (TPR): $TPR")

    TNR = 0
    if N != 0
        TNR = TN / N
    end
    println("specificity, selectivity or true negative rate (TNR): $TNR")

    PPV = 0
    if (TP + FP) != 0
        PPV = TP / (TP + FP)
    end 
    println("precision or positive predictive value (PPV): $PPV")

    NPV = 0
    if (TN + FN) != 0
       NPV = TN / (TN + FN)
    end
    println("negative predictive value (NPV): $NPV")

    FNR = 1 - TPR
    println("miss rate or false negative rate (FNR): $FNR")

    FPR = 1 - TNR
    println("Fall-out or false positive rate (FPR): $FPR")

    FDR = 1 - PPV
    println("false discovery rate: $FDR")

    FOR = 1 - NPV
    println("false omission rate: $FOR")

    ACC = (TP + TN) / total
    println("accuracy: $ACC\n")

    Fone = 0
    if (2 * TP + FP + FN) != 0
        Fone = (2 * TP) / (2 * TP + FP + FN)
    end
    println("F1 Score: $Fone")
   
    MCC = 0
    divideBy = ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))^(1/2)
    if divideBy != 0
        MCC = ((TP * TN) - (FP * FN)) / divideBy
    end
    println("Mathews correlation coefficient (MCC): $MCC")

    BM = TPR + TNR -1
    println("Informedness or Bookmaker Informedness (BM): $BM")

    MK = PPV + NPV -1
    println("Markedness (MK): $MK\n")

    println("Expected 0 predicted 0: $expected_zero_got_zero")
    println("Expected 0 predicted 1: $expected_zero_got_one")
    println("Expected 0 predicted 2: $expected_zero_got_two")
    println("Expected 1 predicted 0: $expected_one_got_zero")
    println("Expected 1 predicted 1: $expected_one_got_one")
    println("Expected 1 predicted 2: $expected_one_got_two")
    println("Expected 2 predicted 0: $expected_two_got_zero")
    println("Expected 2 predicted 1: $expected_two_got_one")
    println("Expected 2 predicted 2: $expected_two_got_two\n\n")


    println("Patient\tNode\tExpected\tActual\n")
    for patientname in keys(errors)
        patient = errors[patientname]
        for genename in keys(patient)
            values = errors[patientname][genename]
            println("$patientname\t$genename\t", values[1], "\t", values[2])
        end
    end
end

end
