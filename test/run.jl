using FactCheck

facts("Interacion Types") do

    tests = ["and", "andneg", "or"]#, "No_Positive_Feedback_Loop_Prototype", "Positive_Feedback_Loop_Prototype" ]

    for (index, value) in enumerate(tests)
        a=readstring(`julia bin/runInference.jl --onenormal ./test/files/networks/$value.tsv ./test/files/observations/$value.tsv ./test/files/results/$value.tsv ./data/db_id_to_name_mapping.txt ./data/key_outputs.tsv -v`)
        std_out=readstring(`diff ./test/files/expected_results/$value.tsv ./test/files/results/$value.tsv`)

        context("$value") do
            @fact std_out --> ""
        end
    end
end
