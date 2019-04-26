module Pi

using DataFrames
using CSV

export readPiFile

struct ModelANDParents
   relation::Any
   parents::Any
end

struct ModelORParents
   relation::Any
   posParents::Any
   negParents::Any
end

struct PIparents
   andPosParents::Array
   andNegParents::Array
   orPosParents::Array
   orNegParents::Array
end

# Assumptions:
## If there is only one parent it will be an AND relation

function readFile(fname)
    df = CSV.read(fname,
                  delim="\t",
                  nullable=false,
                  header=["parentName", "childName", "posneg", "andor"],
                  types=[String, String, Int, Int])

    PIs = Dict{AbstractString,Any}()

    for row in eachrow(df)
        posnegbool = row[:posneg] == 1 ? true : false
        andBool = row[:andor] == 0 ? true : false

        if haskey(PIs, row[:childName])
            if posnegbool
                parents = andBool ?
                             PIs[row[:childName]].andPosParents :
                             PIs[row[:childName]].orPosParents
            else
                parents = andBool ?
                             PIs[row[:childName]].andNegParents :
                             PIs[row[:childName]].orNegParents
            end

            push!(parents, row[:parentName])
        else
            posParents = AbstractString[]
            negParents = AbstractString[]
            if posnegbool
                push!(posParents, row[:parentName])
            else
                push!(negParents, row[:parentName])
            end

            node = PIparents(andBool ? posParents : AbstractString[],
                             andBool ? negParents : AbstractString[],
                             andBool ? AbstractString[] : posParents,
                             andBool ? AbstractString[] : negParents)

            PIs[row[:childName]] = node
        end
    end

    for row in eachrow(df)
        if !haskey(PIs, row[:parentName])
            PIs[row[:parentName]] = PIparents(AbstractString[],AbstractString[],AbstractString[],AbstractString[])
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
                orNodeName = string(nodeName, "_OR")
                nodes[orNodeName] = ModelORParents("OR", PIs[nodeName].orPosParents, PIs[nodeName].orNegParents)
                push!(parents, orNodeName)
            end

            if ((length(PIs[nodeName].andNegParents) > 1) && (length(parents) == 0))
                pseudoParent = "PSEUDONODE_PARENT_$nodeName"
                push!(PIs[nodeName].andPosParents, pseudoParent)
                nodes[pseudoParent] = ModelANDParents("ROOT", AbstractString[])
                parents = PIs[nodeName].andPosParents
            end

            negParents = PIs[nodeName].andNegParents
            if length(parents) == 1 && length(negParents) == 1
                nodes[nodeName] = ModelANDParents("ANDNEG", union(parents, negParents))
            elseif length(parents) == 0 && length(negParents) == 1
                nodes[nodeName] = ModelANDParents("NEG", negParents)
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
                       newChildNodeName = string(childNodeName, "_PSEUDONODE_ ", parenti)
                       nodes[newChildNodeName]= ModelANDParents("AND", AbstractString[childNodeName, parenti])
                       childNodeName = newChildNodeName
                   end
               end

               for i = 1:length(negParents)
                   parenti = negParents[i]
                   if i == length(negParents)
                       nodes[nodeName]= ModelANDParents("AND", AbstractString[childNodeName, parenti])
                   else
                       newChildNodeName = string(childNodeName, "_PSEUDONODE_", parenti)
                       nodes[newChildNodeName]= ModelANDParents("ANDNEG", AbstractString[childNodeName, parenti])
                       childNodeName = newChildNodeName
                   end
               end
            end
        end

    end

    return nodes
end

end
