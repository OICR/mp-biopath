
println(LOAD_PATH)

function analysis(fname)
    println(fname)

    f1= open(fname)
    data = readlines(f1)
    for line in data
        println(line)
    end

    close(f1)
end
println("Running analysis on file:")

if length(ARGS) == 0
    println("Please provide the path to the pairwise interaction file in the first line")
else 
    println("Performing analysis on pairwise interaction file:")
    println(ARGS[1])
 #   analysis(ARGS[1])
end

