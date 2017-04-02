using FactCheck

facts("Interacion Types") do

    @fact 1 --> 1

    context("AND") do
        @fact 1 --> 1
    end

    context("OR") do
        @fact 1 --> not(2)
        @fact 2 --> not(isodd)
    end

    context("AND with Negative interaction") do
        for i=1:5
            @fact i --> i
        end
    end
end
