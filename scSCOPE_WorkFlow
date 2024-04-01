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
expr.data <- fread(filename)
ex.data <- dataFiltering(expr.data)

dir.create("correlation", showWarnings = FALSE)
dir.create("diffCox", showWarnings = FALSE)
dir.create("results_sparse", showWarnings = FALSE)
dir.create("results_sparse/corSparse", showWarnings = FALSE)

correlation_null(ex.data)
diffCorrelation_null(ex.data)
sparse_lasso(ex.data, ITERS_L)
coexpression(ex.data)
pathwayEnrichment(pathgeneId = "ensembl_gene_id", organism = "mmusculus", pathwayDatabase = "KEGG")
markergenes(ex.data, FC_thres = 0.5)
