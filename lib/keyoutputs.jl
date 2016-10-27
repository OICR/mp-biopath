module Keyoutputs

function getNodes()
    data = readlines("./data/keyoutputs.tsv")
    shift!(data)

    keyoutputs = Dict{AbstractString,Any}()

    for line in data
         linefields = split(chomp(line), '\t')
         keyoutputs[linefields[4]] = 1
    end

    return Set(collect(keys(keyoutputs)))
end

end
