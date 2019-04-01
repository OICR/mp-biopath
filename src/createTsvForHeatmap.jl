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
       mergedResultsPathwayLevelFile = "$(resultsFolder)pathway-level-results.tsv"

       outfile = open(mergedResultsFile, "w")
       outPathwayFile = open(mergedResultsPathwayLevelFile, "w")

       if verbose
           println("creating: $mergedResultsFile and $mergedResultsPathwayLevelFile")
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

                  write(outPathwayFile, "nodes\t")
                  for col in columns
                      if col != :node
                         write(outPathwayFile, "$col\t")
                      end
                  end
                  write(outPathwayFile, "pathway_label\n")
               end

               fivePercent = 0.05 * size(columns,1)
               
               averageValues = []
               for row in eachrow(df)
                  nodeName = row[:node]
                  if haskey(keyoutputs, pathwayName) && haskey(keyoutputs[pathwayName], nodeName)
                      nodeId = keyoutputs[pathwayName][nodeName]

                      countChanged = 0
                      for col in columns
                          if col != :node
                             if row[col] != 1
                                 countChanged += 1 
                             end
                          end
                      end
                      
                      if countChanged >= fivePercent
                          for col in columns
                              if col == :node
                                  write(outfile, "$(pathwayId)-$(nodeId)\t")
                              else
                                  value = row[col]
                                  write(outfile, "$value\t")
                              end
                          end
                          write(outfile, "$shortName\n")
                          push!(averageValues, row)
                      end
                  end
               end

               numElements = length(averageValues)
               if numElements > 0
                   write(outPathwayFile, "$(pathwayId)\t")
                   for col in columns
                      if col != :node
                          sum = 0
                          for row in averageValues
                              sum = sum + log(10, row[Symbol(col)])
                          end
                          average = 10 ^ (sum / numElements)
                          write(outPathwayFile, "$average\t")
                      end
                   end
                   write(outPathwayFile, "$shortName\n")
               end
           end
        end
    end
end

end
