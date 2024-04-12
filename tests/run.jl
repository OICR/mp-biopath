using Test

@testset "Interaction Types" begin

    tests = ["and", "andneg", "or", "noPosLoop", "posLoop", "negLoop"]

    for (index, value) in enumerate(tests)
        a=read(`julia bin/mp-biopath inference --config=./tests/conf/$value.yaml -v`)
        std_out = String(read(`diff ./tests/files/expected_results/$value.tsv ./tests/files/outputs/testset/pathways/$value/results.tsv`))
        @test std_out == ""
    end
end;

#=
facts("Loops") do

    loopTests = ["No_Positive_Feedback_Loop_Prototype", "Positive_Feedback_Loop_Prototype" ]

    for (index, value) in enumerate(loopTests)
        a=read(`julia bin/mp-biopath --copynumber ./test/files/networks/$value.tsv ./test/files/observations/$value.tsv ./test/files/results/$value.tsv ./test/files/db_id_to_name_mapping_loop_tests.txt ./data/key_outputs.tsv -v`)
        std_out=read(`diff ./test/files/expected_results/$value.tsv ./test/files/results/$value/$value.tsv`)

        context("$value") do
            @fact std_out --> ""
        end
    end

end
=#
