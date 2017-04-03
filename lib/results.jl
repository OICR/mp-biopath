module Results

include("valuetostate.jl")

function createcsv(nodesampleresults, columns, resultfilename)
    outfile = open(resultfilename, "w")

    write(outfile, join(columns, "\t"))
    write(outfile, "\n")

    for node in keys(nodesampleresults)
        write(outfile, "$node\t")

        index = 0
        lengthColumns = length(columns)
        for column in columns
            index += 1
            if column == "gene"
                continue
            end
            value = round(nodesampleresults[node][column], 2)
            write(outfile, "$value")
            if lengthColumns == index
                write(outfile, "\n")
            else
                write(outfile, "\t")
            end
        end
    end

    close(outfile)
end

function getResults(resultfilename, downregulatedcutoff, upregulatedcutoff, pgmlab)
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
                state = pgmlab? lineparts[i]: ValueToState.getStateNumber(lineparts[i], downregulatedcutoff, upregulatedcutoff)
                counts[state] = counts[state] + 1
                samplenodestate[column][node] = state
            end
        end
    end

    return Dict("counts" => counts, "samplenodestate" => samplenodestate)
end

function getExpected(expectedfile)
    expected_data = readlines(expectedfile)

    header = shift!(expected_data)
    headerparts = split(chomp(header), "\t")

    samplenodestate = Dict()

    i = 0
    for column in headerparts
        i = i + 1
        if i > 1
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
                node = lineparts[1]
            else
                statename = lineparts[i]
                state = statename == "UP" ? "3":
                            statename == "DOWN"? "1":
                                statename == "NC"? "2": statename
                statename = "fixed"
                counts[state] = counts[state] + 1
                samplenodestate[column][node] = state
            end
        end
    end

    return Dict("counts" => counts, "samplenodestate" => samplenodestate)
end

end
