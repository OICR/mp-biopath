module Pi

export readPiFile

type Node
   posParents::Any
   negParents::Any
   andOr::AbstractString #true is and false is or
end

function readFile(fname)
    f1 = open(fname)
    data = readlines(f1)
    shift!(data)
    shift!(data)

    nodes = Dict{AbstractString,Any}()
   
    for line in data
        (parentName, childName, posneg, andor) = split(chomp(line), '\t')

        posnegbool = posneg == "1"? true: false

        if haskey(nodes, childName)
            parents = posnegbool?
                nodes[childName].posParents:
                nodes[childName].negParents
            push!(parents, parentName)
        else
            posParents = Set{AbstractString}()
            negParents = Set{AbstractString}()
            if posnegbool
                push!(posParents, parentName)
            else
                push!(negParents, parentName)
            end
            node = Node(posParents,
                        negParents,
                        andor == "1"? "AND": "OR")
            nodes[childName] = node
        end
    end

    close(f1)

    return nodes
end

end
