#!/usr/bin/env Rscript

#################################################################################

# Created by: Cadia Chan (Lincoln lab rotation student 2017)
# Date: November 6, 2017

# Script is called from "make_output_files.py" to make visual outputs by
# taking arguments specified in .py file.

# By default the heatmap is always generated.
# Depending on parameters specified by user, the t-SNE plot and Kaplan Meier
# survival curves may be generated.

# Outputs files as .svg files.

#################################################################################
cranMirror <- "https://cloud.r-project.org/"
if (!require(gplots)) {install.packages("gplots", repos=cranMirror)}
if (!require(Rcpp)) {install.packages("Rcpp", repos=cranMirror)}
if (!require(ggplot2)) {install.packages("ggplot2", repos=cranMirror)}
if (!require(dendextend)) {install.packages("dendextendRpp", repos=cranMirror)}

args = commandArgs(trailingOnly = TRUE)

data_file <- args[1]
colour_column_string <- args[2]
output_name <- args[3]
list_pathways <- args[4]
list_colours <- args[5]
visual <- args[6]

library(dendextend)
library(gplots)
# library(ggplot2)

#loads in data file containing donors and nodes/observations
data <- read.table(data_file, header = TRUE, sep = "\t", row.names = 1)
df <- data[1:nrow(data),1:(ncol(data)-1)]
y <- as.matrix(df)
y <- y[] + 0.01

my_palette <- colorRampPalette(c("blue", "white", "red"))(n=199)
row_annotation <- unlist(strsplit(colour_column_string, ", "))

if (grepl("h", visual)) {
    output_htmap <- "gecco_heatmap.svg"
    par(lend = 1)
    svg(output_htmap)
}

pathway_string <- unlist(strsplit(list_pathways, ", "))
colour_string <- unlist(strsplit(list_colours, ", "))

nr = dim(y)[1]
nc = dim(y)[2]

print("creating heatmap")

htmap <- heatmap.2(log10(y),
        scale="none",
        key=TRUE,
        key.title="Pathways",
        keysize =1.00,
        key.xlab="Log Base 10",
        xlab="Donors",
        ylab="Keyoutputs",
        labCol=FALSE,
        symkey=FALSE,
        density.info="none",
        col=my_palette,
        trace="none",
        cexRow=0.01 + 0.1/log10(nr),
        cexCol=0.01 + 0.1/log10(nc),
        margins=c(12,15),
        RowSideColors=row_annotation)
print("created heatmap")
    cd <- htmap$colDendrogram
print("DDD")
print(class(cd))
print("printed cd")
#    pdf("mydendrogram.pdf")
#    plot.new()
#    plot(cd)
#    dev.off()
    dd <- cutree(cd, k=[1:10])
#print("eeeee")
  #  print(dd)
#    num_clusters <- 6
#    print("before")
#    ddcut <- cutree(as.hclust(dd),k=[1:6])#num_clusters)
#    print("after")
 #   print(ddcut)
#    print("end")


#library(dynamicTreeCut)
#cd <- htmap$colDendrogram
#dd <- cutreeDynamic(cd )#, k=1:10)

#k3 <- kmeans(na.omit(t(df)), centers = 6, nstart = 25)
#print(k3$cluster)
#exit()

#d1 <- cophenetic(cd)
#print(d1)

if (grepl("h", visual)) {
    legend("topright",      # location of the legend on the heatmap plot
           legend = pathway_string, # category labels
           col = colour_string,
           border=FALSE,
	   bty="n",
	   lwd = 4,
	   y.intersp = 0.7,
	   cex=0.4)
    svg(output_htmap)
    dev.off() 
}

if (grepl("t", visual)) {
    if (!require("Rtsne")) {install.packages("Rtsne", repos='https://cloud.r-project.org')}

    #sets seed so that tsne plot looks the same when you rerun the script
    set.seed(20)

    #gets information from heatmap dendrogram for use in tsne plot
    ddcut <- cutree(htmap$colDendrogram, k=6)
 #   dd <- as.hclust(htmap$colDendrogram)
    num_clusters <- 6
 #   ddcut <- cutree(as.hclust(dd),k=num_clusters) #k = number of clusters
    #manipulate data file for use in the tsne plot
    m <- data[,1:(ncol(data) - 1)] #removes pathway column
    new_data <- rbind(m, ddcut) #adds in subgrouping row (will become a column once matrix is transposed)
    row.names(new_data)[nrow(new_data)] <- "subgroups" #renames subgrouping row

    #makes colours for tsne plot
    colours <- rainbow(num_clusters)
    names(colours) <- unique(ddcut)

    #transposes dataframe
    x <- t(new_data)
    #make 2d tsne plot
    library(Rtsne)

    output_tsne <- "/home/awright/gecco/gecco_tsne.svg" #paste("./", output_name, "_tsne.svg", sep = "")
    par(lend = 1)
    svg(output_tsne)

    tsne_plot <- Rtsne(as.matrix(unique(x[,1:(ncol(x) - 1)])), dims = 2, perplexity =50, max_iter = 5000, check_duplicates =  FALSE)
    plot(tsne_plot$Y, main="tsne", t = 'n')
    text(tsne_plot$Y, labels = x[,ncol(x)], col=colours[x[,ncol(x)] + 1])

    svg(output_tsne)
    dev.off()

    # #make 3d tsne plot
    # library(rgl)

    # output_tsne3d <- paste("./", output_name, "_tsne3d.svg", sep = "")
    # par(lend = 1)
    # svg(output_tsne3d)

    # threed_tsne_plot <- Rtsne(as.matrix(unique(x[,1:113])), perplexity=50, dims=3, check_duplicates = FALSE)
    # plot3d(threed_tsne_plot$Y, labels = x[,114], col=cols[x[,114] + 1])

    # svg(output_tsne3d)
    # dev.off()
}

if (grepl("k", visual)) {

    if (!require('cowplot')) {install.packages("cowplot", repos=cranMirror)}
    if (!require('ggpubr')) {install.packages("ggpubr", repos=cranMirror)}
    if (!require('survival')) {install.packages("survival", repos=cranMirror)}
    if (!require('survminer')) {install.packages("survminer", repos=cranMirror)}

    library(survival)
    # dependencies: cowplot, ggpubr
    library(survminer)

    #reads in toy data set with clinical data
    gecco_test <- read.csv("/home/awright/gecco/gecco_680_kaplan.tsv", sep = "\t", header = F)

    #renames columns names for data frame
    colnames(gecco_test) <- c("donor", "subgroup", "age", "status", "time")
  
    #makes survival object to be plotted and is grouped by SUBGROUP
    gecco.by.subgroup <- survfit(Surv(time, status == 3) ~ subgroup, data = gecco_test, conf.type = "log-log")

    output_kaplan <- "/home/awright/gecco/gecco_kaplan.svg"
    par(lend = 1)
    svg(output_kaplan)

    #makes kaplan meier plots with tables
    surv_plot <- ggsurvplot(gecco.by.subgroup, data = gecco_test, risk.table = TRUE,
               linetype = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), #different types of lines
               #surv.median.line = "v",
               conf.int = FALSE,
               pval = TRUE,
               tables.height = 0.5,
               xlim = c(0, 5), #1825),
               ylim = c(0.8,1),
 #              cumcensor = TRUE,
               #fun = "cumhaz",
               break.time.by = 1
               )
    print(surv_plot)

    svg(output_kaplan)
    dev.off()

    #option to separate each survival curve into their own plot
    # curv_facet <- surv_plot$plot + facet_grid(. ~ subgroup)
}
