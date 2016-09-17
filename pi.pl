module Pi

export readPiFile

type Node
   andPosParents::Any
   andNegParents::Any
   orPosParents::Any
   orNegParents::Any
end

function readFile(fname)
    f1 = open(fname)
    data = readlines(f1)
    shift!(data)
    shift!(data)

    nodes = Dict{AbstractString,Any}()

    #get all nodes that are children first because we know of their and or state
    for line in data
        (parentName, childName, posneg, andor) = split(chomp(line), '\t')

        posnegbool = posneg == "1"? true: false
        andBool = andor == "0"? true: false

        if haskey(nodes, childName)
            if posnegbool 
                parents = andBool?
                             nodes[childName].andPosParents:
                             nodes[childName].orPosParents
            else
                parents = andBool?
                             nodes[childName].andNegParents:
                             nodes[childName].orNegParents
            end

            push!(parents, parentName)
        else
            posParents = Set{AbstractString}()
            negParents = Set{AbstractString}()
            if posnegbool
                push!(posParents, parentName)
            else
                push!(negParents, parentName)
            end
            node = Node(andBool? posParents: Set{AbstractString}(),
                        andBool? negParents: Set{AbstractString}(),
                        andBool? Set{AbstractString}(): posParents,
                        andBool? Set{AbstractString}(): negParents)
                        
            nodes[childName] = node
        end
    end

    #add all nodes that are not children - only parent(s)
    for line in data
        (parentName, childName, posneg, andor) = split(chomp(line), '\t')
        if haskey(nodes, parentName) 
        else
            nodes[parentName] = Node(Set(),Set(),Set(),Set()) 
        end
    end

    close(f1)

    return nodes
end

end
