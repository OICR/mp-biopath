module Pi

export readPiFile

type Node
   posParents::Set{AbstractString}
   negParents::Set{AbstractString}
   andOr::AbstractString #true is and false is or
end

function readFile(fname)
    println(fname)

    f1 = open(fname)
    data = readlines(f1)
    shift!(data)
    shift!(data)

    nodes = Dict{AbstractString,Any}()
   
    for line in data
        (parentName, childName, andor, posneg) = split(chomp(line), '\t')

        posnegbool = posneg == "1"? true: false

        if haskey(nodes, childName)
            parents = posnegbool?
                nodes[childName].posParents:
                nodes[childName].negParents

            println(parents)
 
            push!(parents, parentName)
        else
            node = Node(posnegbool? Set(parentName): Set(),
                        posnegbool? Set(): Set(parentName),
                        andor == "1"? "AND": "OR")
            nodes[childName] = node
        end
    end

    close(f1)

    return nodes
end

end
