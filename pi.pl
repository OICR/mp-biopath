module Pi

export readPiFile

type ModelParents
   relation::Any
   parents::Any
end  

type PIparents
   andPosParents::Any
   andNegParents::Any
   orPosParents::Any
   orNegParents::Any
end

# Assumptions:
## Can't have node with only negative parents
## If there is only one parent it will be an AND relation

function readFile(fname)
    f1 = open(fname)
    data = readlines(f1)
    shift!(data)
    shift!(data)

    PIs = Dict{AbstractString,Any}()

    for line in data
        (parentName, childName, posneg, andor) = split(chomp(line), '\t')

        posnegbool = posneg == "1"? true: false
        andBool = andor == "0"? true: false

        if haskey(PIs, childName)
            if posnegbool 
                parents = andBool?
                             PIs[childName].andPosParents:
                             PIs[childName].orPosParents
            else
                parents = andBool?
                             PIs[childName].andNegParents:
                             PIs[childName].orNegParents
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
            node = PIparents(andBool? posParents: Set{AbstractString}(),
                             andBool? negParents: Set{AbstractString}(),
                             andBool? Set{AbstractString}(): posParents,
                             andBool? Set{AbstractString}(): negParents)
                        
            PIs[childName] = node
        end
    end

    for line in data
        (parentName, childName, posneg, andor) = split(chomp(line), '\t')
        if !haskey(PIs, parentName) 
            PIs[parentName] = PIparents(Set(),Set(),Set(),Set()) 
        end
    end

    nodes = Dict{AbstractString,Any}()
    for nodeName in keys(PIs)
        parents = PIs[nodeName].andPosParents
        if length(PIs[nodeName].orPosParents) > 0 
            if length(PIs[nodeName].andPosParents) == 0 && length(PIs[nodeName].andNegParents) == 0
                nodes[nodeName] = ModelParents("OR", PIs[nodeName].orPosParents)
                continue
            else
                orNodeName = AbstractString(nodeName, "_OR")
                nodes[orNodName] = ModelParents("OR", PIs[nodeName].orPosParents)
                union!(parents, Set{AbstractString}(orNodeName))
            end
        end

        negParents = PIs[nodeName].andNegParents
        if length(parents) + length(negParents) < 3 
           if length(negParents) > 0
                nodeParents = union(parents, negParents)
                nodes[nodeName] = ModelParents("ANDNEG", nodeParents)
           else
                nodes[nodeName] = ModelParents("AND", parents)
           end
        else
           pseudonodeIndex = 1
println("parents\t", parents)

           parentcollection = collect(parents)
println("parentcollection\t", parentcollection) 
           parentone = parentcollection[1]
           parenttwo = parentcollection[2]
println("parents", parentone, "\t", parenttwo)
           pseudonodeName = "$parentone\_$parenttwo\_PSEUDONODE"
println("pseudonodename\t", pseudonodeName)
           nodes[pseudonodeName] = ModelParents("AND", Set(AbstractString(parentone),AbstractString(parenttwo)))
           for i = 3:length(parents)
               newPseudonodeName = "$parentone\_$parenttwo\_PSEUDONODE"
               parenti = parentcollection[i]
               println("prarenti\t", parenti)
               nodes[newPseudoNodeName]= ModelParents("AND", Set{AbstractString}(parenti, pseudoName))
               pseudonodeName = newPseudonodeName
           end 
           nodeIndex = length(parents)
           for i = 1:length(negParents)
               newPseudonodeName = "$parentone\_$parenttwo\_PSEUDONODE"
               nodes[newPseudoNodeName]= ModelParents("ANDNEG", Set{AbstractString}(parentcollection[i], pseudoName))
               pseudonodeName = newPseudonodeName
           end
        end
    end 

    close(f1)

    return nodes
end

end
