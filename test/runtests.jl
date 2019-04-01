using FactCheck

facts("Interaction Types") do

    tests = ["and", "andneg", "or"]

    for (index, value) in enumerate(tests)
        a=readstring(`julia bin/mp-biopath inference --config=./test/conf/$value.conf -v`)
        std_out=readstring(`diff ./test/files/expected_results/$value.tsv ./test/files/results/$value/$value/results.tsv`)

        context("$value") do
            @fact std_out --> ""
        end
    end
end

#=
facts("Loops") do

    loopTests = ["No_Positive_Feedback_Loop_Prototype", "Positive_Feedback_Loop_Prototype" ]

    for (index, value) in enumerate(loopTests)
        a=readstring(`julia bin/mp-biopath --copynumber ./test/files/networks/$value.tsv ./test/files/observations/$value.tsv ./test/files/results/$value.tsv ./test/files/db_id_to_name_mapping_loop_tests.txt ./data/key_outputs.tsv -v`)
        std_out=readstring(`diff ./test/files/expected_results/$value.tsv ./test/files/results/$value/$value.tsv`)

        context("$value") do
            @fact std_out --> ""
        end
    end

end
=#
