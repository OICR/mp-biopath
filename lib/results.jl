module Results

include("valuetostate.jl")

function createcsv(nodesampleresults, columns, resultfilename)
    #resultfilename = join([observationfile, "results"], ".")

    outfile = open(resultfilename, "w")

    write(outfile, join(columns, "\t"))
    write(outfile, "\n")

    for node in keys(nodesampleresults)
        write(outfile, "$node\t")
        for column in columns
            if column == "gene"
                continue
            end
            value = round(nodesampleresults[node][column], 2)
            write(outfile, "$value\t")
        end
        write(outfile, "\n")
    end

    close(outfile)
end

function getResults(observationfile, downregulatedcutoff, upregulatedcutoff)
    resultfilename = join([observationfile, "results"], ".")

    result_data = readlines(resultfilename)

    header = shift!(result_data)
    headerparts = split(chomp(header), "\t")

    samplenodestate = Dict()

    i = 0
    for column in headerparts
        i = i + 1
        if i != 1
            samplenodestate[column] = Dict()
        end
    end

    counts = Dict("1" => 0, "2" => 0, "3" => 0)

    node = ""
    for line in result_data
        lineparts = split(chomp(line), "\t")
        i = 0
        for column in headerparts
            i = i + 1
            if i == 1
                node = lineparts[1]
            else
                state = ValueToState.getStateNumber(lineparts[i], downregulatedcutoff, upregulatedcutoff)
                counts[state] = counts[state] + 1
                samplenodestate[column][node] = state
            end
        end
    end

    return Dict("counts" => counts, "samplenodestate" => samplenodestate)
e

end

function getExpected(expectedfile)
    expected_data = readlines(expectedfile)

    header = shift!(expected_data)
    headerparts = split(chomp(header), "\t")

    samplenodestate = Dict()

    i = 0
    for column in headerparts
        i = i + 1
        if i > 2
            samplenodestate[column] = Dict()
        end
    end

    counts = Dict("1" => 0, "2" => 0, "3" => 0)

    node = ""
    for line in expected_data
        lineparts = split(chomp(line), "\t")
        i = 0
        for column in headerparts
            i = i + 1
            if i == 1
            elseif i == 2
                node = lineparts[2]
            else
                state = lineparts[i]
                counts[state] = counts[state] + 1
                samplenodestate[column][node] = state
            end
        end
    end

    return Dict("counts" => counts, "samplenodestate" => samplenodestate)
end

end
