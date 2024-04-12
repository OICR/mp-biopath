module Expression

using CSV
using DataFrames

"""
# Arguement:
- The filename of the expression file

# Returns
A dataframe of expression data from the file
""" 
function getExpression(filename)
    headerCount = open(filename) do file
        lineCounter = 1
        for line in eachline(file)
             if string(line[1]) == "#"
                lineCounter += 1
             else
                break
             end
        end

        lineCounter
    end

    return CSV.read(filename, delim="\t", header=headerCount)
end


"""
# Arguements:
- expression dataframe
- tissue type

# Returnes
- values for the tissues
"""
function getTissue(df, tissue)
    # Where the Gene Name is the HUGO name
    tissueValues = Dict()
    for row in eachrow(df)
        value = row[tissue]
        if isnothing(value) == false
            tissueValues[row[Symbol("Gene Name")]] = value
        end
    end

    return tissueValues
end

end
