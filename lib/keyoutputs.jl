module Keyoutputs

function getNodes(keyoutputsfile)
    data = readlines(keyoutputsfile)
    shift!(data)

    keyoutputs = Dict{AbstractString,Any}()

    for line in data
         linefields = split(chomp(line), '\t')
         keyoutputs[linefields[5]] = 1
    end

    return Set(collect(keys(keyoutputs)))
end

end
