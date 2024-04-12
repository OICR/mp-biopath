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
		  header=["parent_id",
			  "parent_reactome_id",
			  "parent_name",
			  "parent_type",
			  "child_id",
			  "child_reactome_id",
			  "child_name",
			  "child_type",
			  "polarity",
			  "type"]
                  types=[String, String, Int, Int])

    PIs = Dict{AbstractString,Any}()

    for row in eachrow(df)
        posnegbool = row[:polarity] == "POS" ? true : false
        andBool = row[:type] == "AND" ? true : false
        if haskey(PIs, row[:child_id])
            if posnegbool
                parents = andBool ?
                             PIs[row[:child_id]].andPosParents :
                             PIs[row[:child_id]].orPosParents
            else
                parents = andBool ?
                             PIs[row[:child_id]].andNegParents :
                             PIs[row[:child_id]].orNegParents
            end
            push!(parents, row[:parent_id])
        else
            posParents = AbstractString[]
            negParents = AbstractString[]
            if posnegbool
                push!(posParents, row[:parent_id])
            else
                push!(negParents, row[:parent_id])
            end

            node = PIparents(andBool ? posParents : AbstractString[],
                             andBool ? negParents : AbstractString[],
                             andBool ? AbstractString[] : posParents,
                             andBool ? AbstractString[] : negParents)
            PIs[row[:child_id]] = node
        end
    end

    for row in eachrow(df)
        if !haskey(PIs, row[:parent_id])
            PIs[row[:parent_id]] = PIparents(AbstractString[],AbstractString[],AbstractString[],AbstractString[])
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
                       newChildNodeName = string(childNodeName, "_PSEUDONODE_", parenti)
                       nodes[newChildNodeName]= ModelANDParents("AND", AbstractString[childNodeName, parenti])
                       childNodeName = newChildNodeName
                   end
               end

               for i = 1:length(negParents)
                   parenti = negParents[i]
                   if i == length(negParents)
                       nodes[nodeName]= ModelANDParents("ANDNEG", AbstractString[childNodeName, parenti])
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

function get_nodes(pi_network)
    nodes = Set{String}()
    for node in keys(pi_network)
        if !occursin("PSEUDONODE", node)
            push!(nodes, node)
        end
    end
    return nodes
end

function get_num_edges(pi_network)
    number_of_edges = 0
    for node in keys(pi_network)
        node_obj = pi_network[node]
        node_type = typeof(node_obj)
        for field in fieldnames(node_type)
           if field == :parents
               for node in node_obj.parents
                   if !occursin("PSEUDONODE", node)
                       number_of_edges += 1
                   end                   
               end
           end
           if field == :posParents
               for node in node_obj.posParents
                   if !occursin("PSEUDONODE", node)
                       number_of_edges += 1
                   end                   
               end
           end
           if field == :negParents
               for node in node_obj.negParents
                   if !occursin("PSEUDONODE", node)
                       number_of_edges += 1
                   end                   
               end
           end
        end
    end
    return number_of_edges
end
end
