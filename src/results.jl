module Results

using Printf
include("valuetostate.jl")

function createcsv(nodeSampleResults, columns, resultfilename)
    outfile = open(resultfilename, "w")

    sortedColumns = sort(collect(columns))
    sortedColumnsFull = copy(sortedColumns)
    pushfirst!(sortedColumnsFull, "node")

    write(outfile, join(sortedColumnsFull, "\t"))
    write(outfile, "\n")

    for node in keys(nodeSampleResults)
        if occursin("PSEUDONODE", node)
            continue
        end
        write(outfile, "$node\t")

        index = 0
        lengthColumns = length(columns)
        for column in sortedColumns
            index += 1
            value = @sprintf "%.2f" nodeSampleResults[node][column]
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

function outputAllResults(nodesampleresults, columns, resultfilename)
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
            (x, x_bar) = nodesampleresults[node][column]
            write(outfile, "$x:$x_bar")
            if lengthColumns == index
                write(outfile, "\n")
            else
                write(outfile, "\t")
            end
        end
    end

    close(outfile)
end


function getResults(resultfilename)
    result_data = readlines(resultfilename)

    header = popfirst!(result_data)
    headerparts = split(chomp(header), "\t")

    samplenodestate = Dict()

    i = 0
    for column in headerparts
        i = i + 1
        if i != 1
            samplenodestate[column] = Dict()
        end
    end

    counts = Dict("0" => 0, "1" => 0, "2" => 0)

    node = ""
    for line in result_data
        lineparts = split(chomp(line), "\t")
        i = 0
        for column in headerparts
            i = i + 1
            if i == 1
                node = lineparts[1]
            else
                samplenodestate[column][node] = lineparts[i]
            end
        end
    end

    return samplenodestate
end

function getExpected(expectedfile)
    expected_data = readlines(expectedfile)

    header = popfirst!(expected_data)
    headerparts = split(chomp(header), "\t")

    samplenodestate = Dict()

    i = 0
    for column in headerparts
        i = i + 1
        if i > 1
            samplenodestate[column] = Dict()
        end
    end

    counts = Dict("0" => 0, "1" => 0, "2" => 0)

    node = ""
    for line in expected_data
        lineparts = split(chomp(line), "\t")
        i = 0
        for column in headerparts
            i = i + 1
            if i == 1
                node = lineparts[1]
            else
                state = lineparts[i]
                if state != "-999"
                    counts[state] = counts[state] + 1
                    samplenodestate[column][node] = state
                end
            end
        end
    end

    return Dict("counts" => counts, "samplenodestate" => samplenodestate)
end

end
