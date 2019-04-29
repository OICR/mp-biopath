library(ggplot2)
library(RColorBrewer)
library(d3heatmap)
library(htmlwidgets)
library(pheatmap)

#reads in the synthetic lethal output file
data_file <- read.csv("/home/c2chan/geneList-2keyoutput.tsv", sep = "\t", header = TRUE)
mutated_genes <- read.csv("/home/c2chan/nodes_in_NEW_ROCK.txt", sep = "\t", header = FALSE)

list_genes <- as.vector(mutated_genes$V2)

#function for calculating the difference between synthetic lethal effect and the 
#individual node effect that is greatest between the pair
# calcDiff <- function(node_one, node_two, syn_let) {
#   if(node_one < node_two) {
#     compare <- node_one
#   }else compare <- node_two
#   diff <- 1-(syn_let/compare)
#   return (diff)
# }
calcDiff <- function(node_one, node_two, syn_eff) {
  diff <- abs(node_one - node_two)
  return(syn_eff + diff)
}
#calculates the difference between synthetic lethal effect and the 
#individual node effect that is greatest between the pair
data_file$synergy_value <- mapply(calcDiff, data_file$node_one_value, 
                        data_file$node_two_value, 
                        data_file$effected_node_value)

output_nodes <- unique(as.vector(data_file$effected_node))
#13/130 in mutated genes are root nodes in Rho/ROCK pathway


#creates a dataframe and heat map for each unique effected node (output_node)
for(i in 1:length(output_nodes)) {

   df <- data_file[data_file$effected_node == output_nodes[i],]
   df.col <- as.vector(unique(df$node_one))
   df.row <- list_genes
   new.df <- data.frame(matrix(ncol = length(df.col), nrow = length(df.row)))
   rownames(new.df) <- list_genes
   colnames(new.df) <- df.col
   for (gene.x in list_genes) {
     if (gene.x %in% df$node_one) {
       for (gene.y in df.col) {
         if (gene.y %in% df$node_two) {
           val <- df[df$node_one == gene.x & df$node_two == gene.y, "synergy_value"]
           index.x <- which(rownames(new.df) == gene.x) 
           index.y <- which(names(new.df) == gene.y)
           if(length(val) != 0) { new.df[index.x, index.y] <- val }
       }
     }
       } else if (gene.x %in% df$node_two) {
       for (gene.y in df.col) {
         if (gene.y %in% df$node_one) {
           val <- df[df$node_two == gene.x & df$node_one == gene.y, "synergy_value"]
           index.x <- which(rownames(new.df) == gene.x)
           index.y <- which(names(new.df) == gene.y)
           if(length(val) != 0) { new.df[index.x, index.y] <- val }
         }
       }
     }
   }
   if (output_nodes[i] == "5668934") {
     output_nodes[i] <- "p-T19-MRLC"
   } 
   if (output_nodes[i] == "419195") {
     output_nodes[i] <- "p-T19,S20-MRLC"
   }
   htmap <- d3heatmap(new.df, Rowv = FALSE, Colv = FALSE, dendrogram = "none", 
                      colors = rev("YlGnBu"), digits = 4, na.rm = TRUE,
                      xaxis_height = 300, yaxis_width = 300,
                      yaxis_font_size = 12, labRow = as.vector(mutated_genes$V1))
   print (output_nodes[i])
   # ?saveWidget
   path <- "/home/c2chan/"
   file_name <- output_nodes[i]
   phtmap <- pheatmap(new.df, colorRampPalette(rev(brewer.pal(n = 5, name = "YlGnBu")))(100),
            cluster_rows = FALSE, cluster_cols = FALSE, border_color = "grey100",
            cellwidth = 5,
            fontsize_col = 5, labels_row = as.vector(mutated_genes$V1), main = output_nodes[i],
            filename = paste0(path, file_name, ".pdf", sep = ""))
   saveWidget(htmap, paste0(path, file_name, ".html", sep = ""))
   # assign(output_nodes[i], htmap)
}
# ?pheatmap
# pheatmap(new.df, colorRampPalette(rev(brewer.pal(n = 5, name = "YlGnBu")))(100),
#          cluster_rows = FALSE, cluster_cols = FALSE, border_color = "grey100",
#          cellwidth = 5,
#          fontsize_col = 5, labels_row = as.vector(mutated_genes$V1),
#          filename = "/home/c2chan/test.pdf")
# for (gene.y in df.col) {
#   if (gene.y %in% df$node_one) {
#     for (gene.x in df.row) {
#       if (gene.x %in% df$node_two) {
#         val <- df[df$node_two == gene.x & df$node_one == gene.y, "synergy_value"]
#         index.x <- which(names(new.df) == gene.x) 
#         index.y <- which(names(new.df) == gene.y)
#         if(length(val) != 0) {
#          new.df[index.x, index.y] <- val
#         } else {new.df[gene.y,gene.x] <- 0} 
#       }
#       else { new.df[gene.y,gene.x] <- 0}
#     }
#   }
# }
# ?scale_color_brewer
# MYLK <- ggplot(data = df, aes(x = df$node_one, y = df$node_two, fill = df$diff)) +
#   geom_tile() + scale_fill_distiller(type = "seq", palette = "Blues", direction = 1)
# MYLK
# 
# gene_pairs <- unique(as.data.frame(with(data_file, paste(node_one,",", node_two))))
# ?d3heatmap
# MYLK <- d3heatmap(`p-S1208_S1759-MYLK(1-1914)_[cytosol]_405886`, Rowv = FALSE, Colv = FALSE, dendrogram = "none", 
#           colors = "Blues", symm = TRUE, digits = 4)
# MYLK

