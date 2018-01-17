module CreateTsvForHeatmap

using CSV
using DataFrames

include("idMap.jl")
include("keyoutputs.jl")
include("pathwayList.jl")

function run(resultsFolder, pathwayListFile, keyOutputsFile, verbose)
    keyoutputs = KeyOutputs.getKeyoutputs(keyOutputsFile)
    pathwayList = PathwayList.getPathwayList(pathwayListFile)

    printHeader = false
    if !isdir(resultsFolder)
       println("$resultsFolder is not a directory")
    else
       mergedResultsFile = "$(resultsFolder)results.tsv"
       outfile = open(mergedResultsFile, "w")

       if verbose
           println("creating: $mergedResultsFile")
       end

       dirContents = readdir(resultsFolder)
       columns = ()
       for pathwayName in dirContents
           path = "$resultsFolder$pathwayName"
           if verbose
               println("Processing: $path")
           end 

           if isdir(path)
               pathwayId = get(pathwayList[pathwayName][:pathway_id])
               colour = get(pathwayList[pathwayName][:colour])
               shortName = get(pathwayList[pathwayName][:short_name])
               resultsFile = "$path/results.tsv"
               df = CSV.read(resultsFile,
                             delim="\t",
                             datarow=2,
                             quotechar="\\",
                             nullable=false)

               if printHeader == false
                  printHeader = true
                  columns = names(df)
                  write(outfile, "nodes\t")
                  for col in columns
                      if col != :node
                         write(outfile, "$col\t")
                      end
                  end
                  write(outfile, "pathway_label\n")
               end

               for row in eachrow(df)
                  nodeName = row[:node]
                  if haskey(keyoutputs, pathwayName) && haskey(keyoutputs[pathwayName], nodeName)
                      nodeId = keyoutputs[pathwayName][nodeName]

                      found = false
                      currentValue = -1
                      for col in columns
                          if col != :node
                              if currentValue == -1
                                  currentValue = row[col]
                              end 
                              if row[col] != currentValue
                                  value = row[col]
                                  found = true
                                  break
                              end
                          end
                      end

                      if found == true
                          for col in columns
                              if col == :node
                                  write(outfile, "$(pathwayId)-$(nodeId)\t")
                              else
                                  value = row[col]
                                  write(outfile, "$value\t")
                              end
                          end
                          write(outfile, "$shortName\n")
                      end
                  end
               end
           end
        end
    end
end

end
