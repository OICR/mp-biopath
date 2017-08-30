module CoExpress

include("dbidnamemapping.jl")

using Requests

import Requests: get

function getForPiNodes(pinodes, dbidfile)

    #This returns Hugo Gene ids (what we mainly go by)
    nodetogene = DbIdNameMapping.nodeToGene(dbidfile)

    #responses use entrez gene id
    hugotoentrezid = Dict()
    entrezidtohugo = Dict()
    genes = Dict()
    for node in pinodes
        if haskey(nodetogene, node) == true
            gene = nodetogene[node]
            if haskey(genes, gene) == false
               #gets highly coexpressed genes
               url = "http://coxpresdb.jp/cgi-bin/api2.cgi?gene=$gene&type=mr&cutoff=20&db=hsa2.v12-08"
               ce = get(url)

               sc = statuscode(ce)
               if sc != 200
                   println("got status code $sc when querying $url")
               else
                   response = Requests.json(ce)
                   request = response["request"]
                   entrezid = request["entrez_gene_id"]
                   hugotoentrezid[gene] = entrezid
                   entrezidtohugo[entrezid] = gene
                   genes[gene] = response["results"]
               end
            end
        end
    end


    println(hugotoentrezid)

    genetogene = Dict()
    for (gene, coexpress) in genes
        println(gene)
        for co in coexpress
            cogene = co["gene"]
            if haskey(entrezidtohugo, cogene) == true
                entrezid = entrezidtohugo[cogene]
                println("\t$entrezid")
            end
        end
    end
    exit()

    return [genestoentrezid, genes, genetogene]
end

end
