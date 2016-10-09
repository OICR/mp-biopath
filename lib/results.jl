module Results

function createcsv(nodesampleresults, columns, observationfile)
    resultfilename = join([observationfile, "results"], ".")

    outfile = open(resultfilename, "w")

    write(outfile, join(columns, "\t"))
    write(outfile, "\n")

    for node in keys(nodesampleresults)
        write(outfile, "$node\t")
        for column in columns
            if column == "Gene"
                continue
            end
            value = round(nodesampleresults[node][column], 2)
            write(outfile, "$value\t")
        end
        write(outfile, "\n")
    end

    close(outfile)
end

end

