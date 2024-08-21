library(data.table)
library(caret)
library(caTools)
library(glmnetUtils)
library(doParallel)
library(tictoc)
library(stringr)
library(glmnet)
library(Seurat)
library(dplyr)
library(qlcMatrix)
library(data.table)
library(stringr)
library(splitstackshape)
library(WebGestaltR)
library(KEGGREST)
library(Seurat)

registerDoParallel(2)


START_SEED = 2222
filename = "scope.csv"

PREFIX = str_remove(filename, ".csv")
PREFIX = str_remove(PREFIX, "./")

if(!file.exists(filename)){
  stop(paste0("Input file ", filename, " does not exist! Please check the file path and retry."))
}

print(paste0("Loading file: ", filename))
ex.data1 <- fread(filename)
ex.data <- dataFiltering(ex.data1)

dir.create("correlation", showWarnings = FALSE)
dir.create("diffCox", showWarnings = FALSE)
dir.create("results_sparse", showWarnings = FALSE)
dir.create("results_sparse/corSparse", showWarnings = FALSE)
dir.create("results_sparse/corSparse/all_path", showWarnings = FALSE)


#Now you can find marker genes using any of the following functions: 

findMarkerGenes(). #Can find marker genes between two clusters or one cluster vs all other cells
findAllMarkerGenes() #Can find marker genes for each cluster compared with all other cells
findAllMarkerGenes1vs1() #Can find marker genes between each pair of clusters
geneNetworkPlots() #Can create geneNetwork Plots for a marker gene of interest in selected pathways.
pathwayNetworkPlots() #Can create pathwayNetwork Plots for a pathway in the cluster of interst compared with another cluster or all other cells. 
