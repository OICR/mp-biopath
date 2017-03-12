module Pi

export readPiFile

type ModelANDParents
   relation::Any
   parents::Any
end

type ModelORParents
   relation::Any
   posParents::Any
   negParents::Any
end

type PIparents
   andPosParents::Array
   andNegParents::Array
   orPosParents::Array
   orNegParents::Array
end

# Assumptions:
## Can't have node with only negative parents
## If there is only one parent it will be an AND relation

function readFile(fname)
    f1 = open(fname)
    data = readlines(f1)

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
            posParents = AbstractString[]
            negParents = AbstractString[]
            if posnegbool
                push!(posParents, parentName)
            else
                push!(negParents, parentName)
            end

            node = PIparents(andBool? posParents: AbstractString[],
                             andBool? negParents: AbstractString[],
                             andBool? AbstractString[]: posParents,
                             andBool? AbstractString[]: negParents)

            PIs[childName] = node
        end
    end

    for line in data
        (parentName, childName, posneg, andor) = split(chomp(line), '\t')
        if !haskey(PIs, parentName)
            PIs[parentName] = PIparents(AbstractString[],AbstractString[],AbstractString[],AbstractString[])
        end
    end

    nodes = Dict{AbstractString,Any}()
    for nodeName in keys(PIs)
        parents = PIs[nodeName].andPosParents

        if (length(PIs[nodeName].andPosParents) + length(PIs[nodeName].andNegParents)) == 0
            if (length(PIs[nodeName].orPosParents) + length(PIs[nodeName].orNegParents)) > 0
                nodes[nodeName] = ModelORParents("OR", PIs[nodeName].orPosParents, PIs[nodeName].orNegParents)
            else
                nodes[nodeName] = ModelANDParents("ROOT", parents)
            end
        else
            if (length(PIs[nodeName].orPosParents) + length(PIs[nodeName].orNegParents)) > 0
                orNodeName = "$nodeName\_OR"
                nodes[orNodeName] = ModelORParents("OR", PIs[nodeName].orPosParents, PIs[nodeName].orNegParents)
                push!(parents, orNodeName)
            end

            negParents = PIs[nodeName].andNegParents
            if length(parents) == 1 && length(negParents) == 1
                nodes[nodeName] = ModelANDParents("ANDNEG", union(parents, negParents))
            elseif length(parents) == 1 && length(negParents) == 0
                nodes[nodeName] = ModelANDParents("AND", parents)
            else
               pseudonodeIndex = 1
               childNodeName = parents[1]
               for i = 2:length(parents)
                   parenti = parents[i]
                   if length(negParents) == 0 && i == length(parents)
                       nodes[nodeName]= ModelANDParents("AND", AbstractString[childNodeName, parenti])
                       continue
                   else
                       newChildNodeName = "$childNodeName\_PSEUDONODE\_$parenti"
                       nodes[newChildNodeName]= ModelANDParents("AND", AbstractString[childNodeName, parenti])
                       childNodeName = newChildNodeName
                   end
               end

               for i = 1:length(negParents)
                   parenti = negParents[i]
                   if i == length(negParents)
                       nodes[nodeName]= ModelANDParents("AND", AbstractString[childNodeName, parenti])
                   else
                       newChildNodeName = "$childNodeName\_PSEUDONODE\_$parenti"
                       nodes[newChildNodeName]= ModelANDParents("ANDNEG", AbstractString[childNodeName, parenti])
                       childNodeName = newChildNodeName
                   end
               end
            end
        end

    end

    close(f1)

    return nodes
end

end
