module Pi

export readPiFile

type ModelParents
   relation::Any
   parents::Any
end  

type PIparents
   andPosParents::Array
   andNegParents::Array
   orPosParents::Array
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
                parents = PIs[childName].andNegParents
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
                             andBool? AbstractString[]: posParents)
            
            PIs[childName] = node
        end
    end

    for line in data
        (parentName, childName, posneg, andor) = split(chomp(line), '\t')
        if !haskey(PIs, parentName) 
            PIs[parentName] = PIparents(AbstractString[],AbstractString[],AbstractString[]) 
        end
    end

    nodes = Dict{AbstractString,Any}()
    for nodeName in keys(PIs)
        parents = PIs[nodeName].andPosParents

        if(length(PIs[nodeName].andPosParents) + length(PIs[nodeName].andNegParents)) == 0
            if length(PIs[nodeName].orPosParents) > 0
                nodes[nodeName] = ModelParents("OR", PIs[nodeName].orPosParents)
            else
                nodes[nodeName] = ModelParents("ROOT", parents)
            end
        else
            if length(PIs[nodeName].orPosParents) > 0 
                orNodeName = "$nodeName\_OR"
                nodes[orNodeName] = ModelParents("OR", PIs[nodeName].orPosParents)
                push!(parents, orNodeName)
            end
    
            negParents = PIs[nodeName].andNegParents
            if length(parents) == 1 && length(negParents) == 1
                nodes[nodeName] = ModelParents("ANDNEG", union(parents, negParents))
            else
               pseudonodeIndex = 1
               childNodeName = parents[1]
               for i = 2:length(parents)
                   parenti = parents[i]
                   if length(negParents) == 0 && i == length(parents)
                       nodes[nodeName]= ModelParents("AND", AbstractString[childNodeName, parenti])
                       continue
                   else
                       newChildNodeName = "$childNodeName\_PSEUDONODE\_$parenti"               
                       nodes[newChildNodeName]= ModelParents("AND", AbstractString[childNodeName, parenti])
                       childNodeName = newChildNodeName
                   end
               end 
    
               for i = 1:length(negParents)
                   parenti = parents[i]
                   if i == length(negParents)
                       nodes[nodeName]= ModelParents("AND", AbstractString[childNodeName, parenti])
                   else
                       newChildNodeName = "$childNodeName\_PSEUDONODE\_$parenti"
                       nodes[newChildNodeName]= ModelParents("ANDNEG", AbstractString[childNodeName, parenti])
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
