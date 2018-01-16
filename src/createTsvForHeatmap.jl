module CreateTsvForHeatmap

include("idMap.jl")
include("keyoutputs.jl")
include("pathwayList.jl")

function run(resultsFolder, pathwayListFile, dbidFile, keyOutputsFile, verbose)
 #   IDmap = IdMap.getIdMap(dbidFile)
  #  keyoutputs = KeyOutputs.getKeyoutputs(keyOutputsFile)
  #  pathwayList = PathwayList.getPathwayList(pathwayListFile)

    if !isdir(resultsFolder)
       println("$resultsFolder is not a directory")
    else
       dirContents = readdir(resultsFolder)
       for content in dirContents
           path = "$resultsFolder$content"
           if isdir(path)
               #if  result.tsv is file read it in and combine results. 
           end
        end
    end
end

end
