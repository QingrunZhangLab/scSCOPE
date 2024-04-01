# Function to flatten the correlation matrix
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.table(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut]
  )
}

dataFiltering <- function(ex.data){
  print(paste0("Data loaded. Dimensions: ", paste0(dim(ex.data), collapse = ",")))
  phen <- ex.data$phenotype
  ex.data <- subset(ex.data, select = -phenotype)
  # Filtering low variance and constant expression values
  print(paste0("Calculating variances..."))
  vars <- sapply(ex.data, var)
  constants <- names(vars)[vars!=0 | is.na(vars)]
  print(paste0("Identified ", sum(vars==0, na.rm=TRUE), " transcripts with 0 variance. Removing..."))
  ex.data <- ex.data[, ..constants]
  print(paste0("Remaining data dimensions: ", paste0(dim(ex.data), collapse = ",")))
  
  print(paste0("Calculating variances..."))
  vars <- sapply(ex.data, var)
  low_vars <- names(vars)[is.na(vars) | vars >= quantile(vars, 0.25, na.rm = TRUE)]
  print(paste0("Identified ", (ncol(ex.data) - length(low_vars) - 2), " transcripts with lowest 25% variance. Removing..."))
  ex.data <- ex.data[, ..low_vars]
  print(paste0("Remaining data dimensions: ", paste0(dim(ex.data), collapse = ",")))
  
  ex.data$phenotype <- phen
  # Taking a look at available phenotypes
  print("Current phenotypes: ")
  print(table(ex.data$phenotype))
  return(ex.data)
}


correlation_null <- function(ex.data, ITERS = 100, CORR_PERCENTILE_THRESHOLD = 0.975, PROBES_PER_ITER = 1000){
  if(file.exists(paste0("correlation/", PREFIX, "_CorrNull.csv"))){
    return()
  }
  null_dist <- c()
  set.seed(START_SEED)
  for(i in 1:ITERS){
    print(paste0("Current trial ", i))
    rnd_smp <- sample(colnames(ex.data)[2:(ncol(ex.data)-1)], PROBES_PER_ITER)
    cor_mat <- cor(ex.data[, ..rnd_smp])
    corrs <- flattenCorrMatrix(cor_mat)
    if(is.null(null_dist)){
      null_dist <- as.data.table(corrs$cor)
    }else{
      null_dist <- cbind(null_dist, corrs$cor)
    }
  }
  print("Correlation Null Distribution Done. Writing all results.")
  fwrite(null_dist, paste0("./correlation/", PREFIX, "_CorrNull.csv"))
}



diffCorrelation_null <- function(ex.data, ITERS = 100, DIFFCORR_PERCENTILE_THRESHOLD = 0.975, PROBES_PER_ITER = 1000){
  
  for(cls_no in 1:length(unique(ex.data$phenotype))){
    if(file.exists(paste0("diffCox/scope", cls_no, "_DiffCorrNull.csv"))){
      next
    }
    expr.data <- ex.data
    expr.data$phenotype <- ifelse(expr.data$phenotype == cls_no, 0,1)
    print("Current phenotypes after substituting: ")
    print(table(expr.data$phenotype))
    cluster <- expr.data[expr.data$phenotype == 0,]
    other <- expr.data[expr.data$phenotype == 1,]
    
    null_dist <- c()
    
    set.seed(START_SEED)
    for(i in 1:ITERS){
      print(paste0("Current trial ", i))
      
      rnd_smp <- sample(colnames(expr.data)[2:(ncol(expr.data)-1)], PROBES_PER_ITER)
      
      cluster_sel <- cluster[, ..rnd_smp]
      other_sel <- other[, ..rnd_smp]
      
      cluster_cor <- cor(cluster_sel)
      other_cor <- cor(other_sel)
      
      cluster_corflat <- flattenCorrMatrix(cluster_cor)
      other_corflat <- flattenCorrMatrix(other_cor)
      
      merged <- merge(cluster_corflat, other_corflat, by = c("row", "column"))
      merged$diffcor <- abs(merged$cor.x - merged$cor.y)
      
      if(is.null(null_dist)){
        null_dist <- as.data.table(merged$diffcor)
      }else{
        null_dist <- cbind(null_dist, merged$diffcor)
      }
    }
    
    print("Writing all results.")
    fwrite(null_dist, paste0("./diffCox/", PREFIX, cls_no, "_DiffCorrNull.csv"))
  }
}



#Now run one vs all logistic regression for each cluster with respect to all data in the dataset.

dir.create("results_sparse", showWarnings = FALSE)

sparse_lasso <- function(ex.data, ITERS_L = 200, FOLDS = 10){
  START_SEED = 2222
  all_clusters = unique(ex.data$phenotype)
  for(x in all_clusters) {
    cls_no = x
    if(file.exists(paste0("results_sparse/", PREFIX, cls_no, "_Summary.csv" ))){
      next}
    
    expr.data <- ex.data
    expr.data$phenotype <- ifelse(expr.data$phenotype == cls_no, 0,1)
    expr.data <- expr.data[, -c(1)]
    results <- c()
    pred_metrics <- c()
    set.seed(START_SEED)
    
    for(k in 1:ITERS_L){
      
      # Splitting data proportional to phenotype composition
      inds  <- sample.split(expr.data$phenotype, SplitRatio = 0.7)
      phen <- expr.data[inds, ]$phenotype
      tmp <- subset(expr.data, select = -phenotype)
      print(paste0("Currently fitting iteration: ", k))
      tic("Iteration time")
      
      # Fitting a glmnet LASSO model for the training data
      train <- as(as.matrix(tmp[inds, ]), "sparseMatrix")
      cvfit <- glmnetUtils::cv.glmnet(train, phen, alpha = 1,
                                      family = "binomial",
                                      nfolds = FOLDS, intercept = FALSE,
                                      parallel = TRUE)
      
      # Obtaining fitted coefficients
      tmp_coeffs <- coef(cvfit)
      print("Fitted...")
      
      # Keeping track of results from each of the LASSO models fit in each iteration
      if(is.null(results)){
        results <- data.table(probe = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1],
                              coefficient_1 = tmp_coeffs@x)
      }else{
        results <- merge(results, data.table(probe = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1],
                                             coefficient = tmp_coeffs@x),
                         all.x=TRUE, all.y=TRUE,
                         by="probe", suffixes = c("",paste0("_", k)))
      }
      
      print("Predicting...")
      tokeep <- which(sapply(tmp, is.numeric))
      test_data <- as.matrix(tmp[!inds, ..tokeep])
      cvpred <- predict(cvfit, test_data)
      
      toc(quiet = FALSE)
      
      print("Calculate prediction metrics...")
      preds <- ifelse(cvpred < 0.5, 0, 1)
      acc <- mean(preds == expr.data$phenotype[!inds])
      sens <- sensitivity(as.factor(preds), as.factor(expr.data$phenotype[!inds]), positive = "1", negative = "0")
      spef <- specificity(as.factor(preds), as.factor(expr.data$phenotype[!inds]), positive = "1", negative = "0")
      
      pred_metrics <- rbindlist(list(pred_metrics,
                                     data.table(model_num = k,
                                                Accuracy = acc,
                                                Sensitivity = sens,
                                                Specificity = spef)))
    }
    print("Combining LASSO results...")
    cof_sum <- data.table(probe = results$probe,
                          model.count = ITERS_L - rowSums(is.na(results[,2:ncol(results)])),
                          min.val = apply(results[,2:ncol(results)], 1, min, na.rm = TRUE),
                          max.val = apply(results[,2:ncol(results)], 1, max, na.rm = TRUE))
    print("Writing all results.")
    
    fwrite(cof_sum, paste0("results_sparse/", PREFIX, cls_no,  "_Summary.csv"))
    fwrite(results, paste0("results_sparse/", PREFIX, cls_no, "_AllCoefs.csv"))
    fwrite(pred_metrics, paste0("results_sparse/", PREFIX, cls_no, "_PredMetrics.csv"))
  }
}

coexpression <- function(ex.data, CORE_CUTOFF = 160, sub_sample = 0.6, corr_iter_cutoff = 80, iter_cor = 100,
                         CORR_PERCENTILE_THRESHOLD = 0.975, DIFFCORR_PERCENTILE_THRESHOLD = 0.975){
  
  for(x in list.files("results_sparse/", pattern = "_Summary.csv")){ 
    if(file.exists(paste0("results_sparse/corSparse/", str_replace(x, "_Summary.csv", "") , "_AllCoreCorrelations.csv"))){
      next
    }
    
    SUMMARY_FILE = x
    y = str_remove(x, "_Summary.csv")
    ################################################################################################################
    EXPR_FILE = "scope.csv"
    CORR_NULL_FILE = "correlation/scope_CorrNull.csv"
    #################################################################################################################
    cls_no = str_replace(y, PREFIX, "")
    #################################################################################################################
    DIFFCORR_NULL_FILE = paste0("diffCox/", y, "_DiffCorrNull.csv")
    print(c(y,DIFFCORR_NULL_FILE))
    
    print(paste0("Loading gene summary file: ", SUMMARY_FILE))
    gene_summary <- fread(paste0("results_sparse/", SUMMARY_FILE))
    
    print("Filtering for core genes...")
    core_genes <- gene_summary$probe[gene_summary$model.count >= CORE_CUTOFF & gene_summary$probe != "(Intercept)"]
    
    if(length(core_genes) < 5){
      core_genes <- gene_summary[order(gene_summary$model.count, decreasing = TRUE), ][1:5,]$probe
      core_genes <- core_genes[!is.na(core_genes)]
    }
    gene_coords <- as.data.table(core_genes)
    colnames(gene_coords) <- "external_gene_name"
    
    print("Reading null distribution of co-expressions...")
    corr_null_dist <- fread(CORR_NULL_FILE)
    positive_cutoff <- c()
    negative_cutoff <- c()
    print("Calculating median percentile from co-expression data...")
    for(i in 1:ncol(corr_null_dist)){
      
      positive_cutoff <- c(positive_cutoff, quantile(corr_null_dist[as.logical(corr_null_dist[, ..i] > 0), ..i][[1]], CORR_PERCENTILE_THRESHOLD)[[1]])
      negative_cutoff <- c(negative_cutoff, -quantile(-corr_null_dist[as.logical(corr_null_dist[, ..i] < 0), ..i][[1]], CORR_PERCENTILE_THRESHOLD)[[1]])
      
    }
    
    pos_cut <- median(positive_cutoff)
    neg_cut <- median(negative_cutoff)
    
    print(c(pos_cut, neg_cut))
    print("Reading null distribution of differential co-expressions...")
    diffcorr_null_dist <- fread(DIFFCORR_NULL_FILE)
    
    diffcorr_cutoff <- c()
    
    print("Calculating median percentile from differential co-expression data...")
    for(i in 1:ncol(diffcorr_null_dist)){
      diffcorr_cutoff <- c(diffcorr_cutoff, quantile(diffcorr_null_dist[, ..i][[1]], DIFFCORR_PERCENTILE_THRESHOLD, na.rm = TRUE)[[1]])
    }
    
    diffcorr_cut <- median(diffcorr_cutoff)
    
    print("Loading expression data to identify secondary genes...")
    expr.data <- ex.data
    
    # Taking a look at available phenotypes
    print("Current phenotypes: ")
    print(table(expr.data$phenotype))
    cluster0 <- expr.data[expr.data$phenotype == cls_no,]
    other0 <- expr.data[expr.data$phenotype !=cls_no,]
    
    corr_res = c()
    for(iter in 1:iter_cor){
      print("Current phenotypes after substituting: ")
      
      print(table(expr.data$phenotype))
      
      cluster <- stratified(cluster0, "phenotype", sub_sample)
      other <- stratified(other0, "phenotype", sub_sample)	
      expr1 <- stratified(expr.data, "phenotype", sub_sample)
      results <- c()
      
      clsM <- as.matrix(cluster[, !c("V1", "phenotype")])
      otherM <- as.matrix(other[, !c("V1", "phenotype")])
      allM <- as.matrix(expr1[,!c("V1", "phenotype")])
      
      clsM_core <- as.matrix(subset(cluster, select = core_genes))
      otherM_core <- as.matrix(subset(other, select = core_genes))
      allM_core <- as.matrix(subset(expr1, select = core_genes))
      
      all_res <- corSparse(allM_core, allM)
      rownames(all_res) <- colnames(allM_core)
      colnames(all_res) <- colnames(allM)
      
      cls_res <- corSparse(clsM_core, clsM)
      other_res <- corSparse(otherM_core, otherM)
      rownames(cls_res) <- colnames(clsM_core)
      colnames(cls_res) <- colnames(clsM)
      rownames(other_res) <- colnames(otherM_core)
      colnames(other_res) <- colnames(otherM)
      
      results_all <- data.table(`Core Gene` = rownames(all_res), all_res)
      results_all <- melt(results_all, id.vars = "Core Gene", value.name = "Correlation", variable.name = "Secondary Gene")
      
      results_cluster <- data.table(`Core Gene` = rownames(cls_res), cls_res)
      results_cluster <- melt(results_cluster, id.vars = "Core Gene", value.name = "Cluster Correlation", variable.name = "Secondary Gene")
      
      results_other <- data.table(`Core Gene` = rownames(other_res), other_res)
      results_other <- melt(results_other, id.vars = "Core Gene", value.name = "Other Correlation", variable.name = "Secondary Gene")
      
      results <- merge(results_cluster, results_other, by = c("Core Gene", "Secondary Gene"), all = TRUE)	
      results <- merge(results, results_all, by = c("Core Gene", "Secondary Gene"), all = TRUE)
      
      
      results$`Differential Correlation` <- abs(results$`Cluster Correlation` - results$`Other Correlation`)
      results <- results[results$`Core Gene` != results$`Secondary Gene`, ]
      
      filtered_results <- results[results$Correlation > pos_cut | results$Correlation < neg_cut | results$`Differential Correlation` > diffcorr_cut,]
      filtered_results[, paste0("iter", iter)] <- as.integer(rowSums(is.na(filtered_results[, 3:5])) != 3)
      curr_final <- subset(filtered_results, select = c("Core Gene", "Secondary Gene", paste0("iter", iter)))
      
      if(is.null(corr_res)){
        corr_res <- curr_final
      }else{
        corr_res <- merge(corr_res, curr_final, by = c("Core Gene", "Secondary Gene"), all.x=TRUE, all.y=TRUE)
      }	
    }
    
    final_corr <- data.table(`Core Gene` = corr_res$`Core Gene`,
                             `Secondary Gene` = corr_res$`Secondary Gene`,
                             model.count = iter_cor - rowSums(is.na(corr_res[,3:ncol(corr_res)])))
    final_corr <- final_corr[final_corr$model.count >= corr_iter_cutoff, ]
    
    res_f <- c()
    
    
    clsM_core <- as.matrix(subset(cluster0, select = unique(final_corr$`Core Gene`)))
    otherM_core <- as.matrix(subset(other0, select = unique(final_corr$`Core Gene`)))
    allM_core <- as.matrix(subset(expr.data, select = unique(final_corr$`Core Gene`)))
    
    clsM_sec <- as.matrix(subset(cluster0, select = unique(final_corr$`Secondary Gene`)))
    otherM_sec <- as.matrix(subset(other0, select = unique(final_corr$`Secondary Gene`)))
    allM_sec <- as.matrix(subset(expr.data, select = unique(final_corr$`Secondary Gene`)))
    
    all_res <- corSparse(allM_core, allM_sec)
    rownames(all_res) <- colnames(allM_core)
    colnames(all_res) <- colnames(allM_sec)
    
    cls_res <- corSparse(clsM_core, clsM_sec)
    other_res <- corSparse(otherM_core, otherM_sec)
    rownames(cls_res) <- colnames(clsM_core)
    colnames(cls_res) <- colnames(clsM_sec)
    rownames(other_res) <- colnames(otherM_core)
    colnames(other_res) <- colnames(otherM_sec)
    
    results_all <- data.table(`Core Gene` = rownames(all_res), all_res)
    results_all <- melt(results_all, id.vars = "Core Gene", value.name = "Correlation", variable.name = "Secondary Gene")
    
    results_cluster <- data.table(`Core Gene` = rownames(cls_res), cls_res)
    results_cluster <- melt(results_cluster, id.vars = "Core Gene", value.name = "Cluster Correlation", variable.name = "Secondary Gene")
    
    results_other <- data.table(`Core Gene` = rownames(other_res), other_res)
    results_other <- melt(results_other, id.vars = "Core Gene", value.name = "Other Correlation", variable.name = "Secondary Gene")
    
    results <- merge(results_cluster, results_other, by = c("Core Gene", "Secondary Gene"), all = TRUE)	
    results <- merge(results, results_all, by = c("Core Gene", "Secondary Gene"), all = TRUE)
    
    print(c("Nrow_Results = ", nrow(results)))
    
    print("Calculating differential co-expression...")
    results$`Differential Correlation` <- abs(results$`Cluster Correlation` - results$`Other Correlation`)
    results <- results[results$`Core Gene` != results$`Secondary Gene`, ]
    
    # Identifying CGNs using secondary genes that are co-expressed or 
    # differentially co-expressed.
    filtered_results <- results[results$Correlation > pos_cut | results$Correlation < neg_cut | results$`Differential Correlation` > diffcorr_cut,]
    
    
    print("Writing all results.")
    
    fwrite(results, paste0("results_sparse/corSparse/", y, "_AllCoreCorrelations.csv"))
    fwrite(filtered_results, paste0("results_sparse/corSparse/", y, "_SecondaryGenes.csv"))
    print("Writing all results.")
  }
}


pathwayEnrichment <- function(pathwayDatabase = c("KEGG", "GOBP", "GOCC", "GOMF"), pathgeneId, organism){
  pathwayDatabase = c("KEGG", "GOBP", "GOCC", "GOMF")
  pdb <- unique(pathwayDatabase)
  dir.create("results_sparse/corSparse/all_path", showWarnings = FALSE)
  
  for(x in list.files("results_sparse/corSparse/", pattern = "_SecondaryGenes.csv")){
    
    SECONDARY_GENE_FILE = x
    full_data <- fread(paste0("results_sparse/corSparse/", SECONDARY_GENE_FILE))
    
    print("File Read")
    
    if("KEGG" %in% pdb){
      if(exists("final_path")){
        rm(final_path)}
      if(!file.exists(paste0("results_sparse/corSparse/all_path/path/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))){
        
        for(curr_core_gene in unique(full_data$`Core Gene`)){
          
          curr_gene_list <- c(curr_core_gene, full_data$`Secondary Gene`[full_data$`Core Gene` == curr_core_gene])
          path_res <- NULL
          tryCatch(path_res <- as.data.table(WebGestaltR(interestGene = curr_gene_list,
                                                         organism = organism,
                                                         enrichDatabase="pathway_KEGG",
                                                         interestGeneType=pathgeneId, referenceSet = "genome", isOutput = FALSE, 
                                                         hostName = "https://www.webgestalt.org")),
                   error = function(c){
                     print("Error in enrichment. Skipping...")
                     path_res <<- NULL
                   })
          
          if(is.null(path_res) || nrow(path_res) == 0){
            print("No pathways enriched. Skipping...")
          } else {
            
            path_res$`Core Gene` <- curr_core_gene
            path_res$Cluster <- str_remove(x, "_Summary.csv")
            
            if(exists("final_path")){
              final_path <- rbindlist(list(final_path, path_res))
            } else{
              final_path <- path_res
            }
          }
        }
      }
      
    }
    
    if("GOCC" %in%pdb){
      if(exists("final_cc")){
        rm(final_cc)}
      if(!file.exists(paste0("results_sparse/corSparse/all_path/cc/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))){
        
        for(curr_core_gene in unique(full_data$`Core Gene`)){
          
          curr_gene_list <- c(curr_core_gene, full_data$`Secondary Gene`[full_data$`Core Gene` == curr_core_gene])
          path_res <- NULL
          tryCatch(path_res <- as.data.table(WebGestaltR(interestGene = curr_gene_list,
                                                         organism = organism,
                                                         enrichDatabase="geneontology_Cellular_Component",
                                                         interestGeneType=pathgeneId, referenceSet = "genome", isOutput = FALSE, 
                                                         hostName = "https://www.webgestalt.org")),
                   error = function(c){
                     print("Error in enrichment. Skipping...")
                     path_res <<- NULL
                   })
          
          if(is.null(path_res) || nrow(path_res) == 0){
            print("No pathways enriched. Skipping...")
          } else {
            
            path_res$`Core Gene` <- curr_core_gene
            path_res$Cluster <- str_remove(x, "_Summary.csv")
            
            if(exists("final_cc")){
              final_cc <- rbindlist(list(final_cc, path_res))
            } else{
              final_cc <- path_res
            }
          }
        }
      }
      
    }
    
    if("GOMF" %in% pdb){
      if(exists("final_mf")){
        rm(final_mf)}
      if(!file.exists(paste0("results_sparse/corSparse/all_path/mf/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))){
        
        for(curr_core_gene in unique(full_data$`Core Gene`)){
          
          curr_gene_list <- c(curr_core_gene, full_data$`Secondary Gene`[full_data$`Core Gene` == curr_core_gene])
          path_res <- NULL
          tryCatch(path_res <- as.data.table(WebGestaltR(interestGene = curr_gene_list,
                                                         organism = organism,
                                                         enrichDatabase="geneontology_Molecular_Function",
                                                         interestGeneType=pathgeneId, referenceSet = "genome", isOutput = FALSE, 
                                                         hostName = "https://www.webgestalt.org")),
                   error = function(c){
                     print("Error in enrichment. Skipping...")
                     path_res <<- NULL
                   })
          
          if(is.null(path_res) || nrow(path_res) == 0){
            print("No pathways enriched. Skipping...")
          } else {
            
            path_res$`Core Gene` <- curr_core_gene
            path_res$Cluster <- str_remove(x, "_Summary.csv")
            
            if(exists("final_mf")){
              final_mf <- rbindlist(list(final_mf, path_res))
            } else{
              final_mf <- path_res
            }
          }
        }
      }
      
    }
   
    if("GOBP" %in% pdb){
      if(exists("final_bp")){
        rm(final_bp)}
      if(!file.exists(paste0("results_sparse/corSparse/all_path/bp/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))){
        
        for(curr_core_gene in unique(full_data$`Core Gene`)){
          
          curr_gene_list <- c(curr_core_gene, full_data$`Secondary Gene`[full_data$`Core Gene` == curr_core_gene])
          path_res <- NULL
          tryCatch(path_res <- as.data.table(WebGestaltR(interestGene = curr_gene_list,
                                                         organism = organism,
                                                         enrichDatabase="geneontology_Biological_Process",
                                                         interestGeneType=pathgeneId, referenceSet = "genome", isOutput = FALSE, 
                                                         hostName = "https://www.webgestalt.org")),
                   error = function(c){
                     print("Error in enrichment. Skipping...")
                     path_res <<- NULL
                   })
          
          if(is.null(path_res) || nrow(path_res) == 0){
            print("No pathways enriched. Skipping...")
          } else {
            
            path_res$`Core Gene` <- curr_core_gene
            path_res$Cluster <- str_remove(x, "_Summary.csv")
            
            if(exists("final_bp")){
              final_bp <- rbindlist(list(final_bp, path_res))
            } else{
              final_bp <- path_res
            }
          }
        }
      }
      
    }

    if(exists("final_bp")){
      fwrite(final_bp, paste0("results_sparse/corSparse/all_path/bp/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))}
    
    if(exists("final_cc")){
      fwrite(final_cc, paste0("results_sparse/corSparse/all_path/cc/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))}
    
    if(exists("final_mf")){
      fwrite(final_mf, paste0("results_sparse/corSparse/all_path/mf/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))}
    
    if(exists("final_path")){
      fwrite(final_path, paste0("results_sparse/corSparse/all_path/path/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))}
  }

}

markergenes <- function(ex.data,  DIFFCORR_PERCENTILE_THRESHOLD = 0.975, FC_thres = 0.5){
  expr.data <- data.frame(as.matrix(ex.data, rownames = 1))
  phen <- expr.data$phenotype
  expr.data <- subset(expr.data, select = -phenotype)
  all_cell <- CreateSeuratObject(t(expr.data), project = "sig")
  Idents(all_cell) <- phen
  ncls = length(unique(phen))
  all_cell$orig.ident = "all"
  
  for(dirr in list.dirs(path = "./results_sparse/corSparse/all_path", recursive = FALSE)){
    
    dir.create(paste0(dirr, "/final"), showWarnings = FALSE)
    
    all_pathways <- c()
    
    for(curr_file in list.files(dirr, pattern = "_Pathways.csv")){
      print(c("Reading file:", paste0(dirr, "/", curr_file)))
      curr_data <- fread(paste0(dirr, "/", curr_file))
      all_pathways <- rbindlist(list(all_pathways, curr_data))
    }
    
    all_corrs <- c()
    
    for(curr_file in list.files("./results_sparse/corSparse", pattern = "_SecondaryGenes.csv")){
      print(c("Reading file:", curr_file))
      curr_data <- fread(paste0("./results_sparse/corSparse/", curr_file))
      
      curr_Cluster <- str_replace(curr_file, "_SecondaryGenes.csv","")
      curr_Cluster <- str_replace(curr_file, PREFIX, "")
      curr_data$Cluster <- curr_Cluster
      
      print(c("Reading Differential Correlation file"))
      diffcorr_null_dist <- fread(paste0("./diffCox/", curr_Cluster, "_DiffCorrNull.csv"))
      diffcorr_cutoff <- c()
      print("Calculating median percentile from differential co-expression data...")
      for(i in 1:ncol(diffcorr_null_dist)){
        diffcorr_cutoff <- c(diffcorr_cutoff, quantile(diffcorr_null_dist[, ..i][[1]], DIFFCORR_PERCENTILE_THRESHOLD, na.rm = TRUE)[[1]])
      }
      
      diffcorr_cut <- median(diffcorr_cutoff)
      
      curr_data$`Differential Correlation Sig` <- diffcorr_cut
      curr_data$`Is Differential Correlation Sig` <- ifelse(curr_data$`Differential Correlation` > curr_data$`Differential Correlation Sig`, "Yes", "No")
      
      all_corrs <- rbindlist(list(all_corrs, curr_data))
    }
    
    final_data <- c()
    for(i in 1:nrow(all_pathways)){
      genes <- unique(str_split(all_pathways$userId[i], ";")[[1]])
      curr_data <- data.table(Cluster = rep(all_pathways$Cluster[i], length(genes)),
                              `KEGG GeneSet` = rep(all_pathways$geneSet[i], length(genes)),
                              `Pathway Name` = rep(all_pathways$description[i], length(genes)),
                              `Core Gene` = rep(all_pathways$`Core Gene`[i], length(genes)),
                              `Secondary Gene` = genes)
      
      final_data <- rbindlist(list(final_data, curr_data))
    }
    
    final_data$`Core Gene` <- all_corrs$`Core Gene`[match(final_data$`Core Gene`, all_corrs$`Core Gene`)]
    final_data$`Secondary Gene` <- all_corrs$`Secondary Gene`[match(final_data$`Secondary Gene`, all_corrs$`Secondary Gene`)]
    final_data$Cluster <- str_remove(final_data$Cluster, "_SecondaryGenes.csv")
    
    final_data <- merge(final_data, all_corrs, by = c("Cluster", "Core Gene", "Secondary Gene"), no.dups = TRUE)
    
    nclss <- unique(final_data$Cluster)
    gl_res <- c()
    for(cls in nclss){
      subb <- final_data[final_data$Cluster == cls, ]
      diff_sig <- subb$'Differential Correlation Sig'[1]
      #subb <- subb[subb$`Is Differential Correlation Sig` == "Yes" | (is.na(subb$`Other Correlation`) & subb$`Cluster Correlation` >= diff_sig) | (is.na(subb$`Other Correlation`) & subb$`Cluster Correlation` <= -diff_sig) | (is.na(subb$`Cluster Correlation`) & subb$`Other Correlation` >=diff_sig)| (is.na(subb$`Cluster Correlation`) & subb$`Other Correlation` <= -diff_sig) | subb$Correlation >= diff_sig | subb$Correlation < -diff_sig, ]
      
      cls = str_remove(cls, "scope")
      cls = as.integer(str_split(cls, "_")[[1]][1])
      c_all <- unique(c(subb$`Core Gene`, subb$`Secondary Gene`))
      ##avg_cls <- AverageExpression(all_cell, features = c_all)
      #avg_all <- AverageExpression(all_cell, features = c_all, group.by = "orig.ident")
      #fold <- FoldChange(all_cell, ident.1 = as.integer(cls), features = c_all)
      
      core <- subb$`Core Gene`
      sec <- subb$`Secondary Gene`
      
      #subb$Core_avgExpr_clus <- avg_cls$RNA[core, cls]
      #subb$Sec_avgExpr_clus <- avg_cls$RNA[sec, cls]
      #subb$Core_avgExpr_all <- avg_all$RNA[core, ]
      #subb$Sec_avgExpr_all <- avg_all$RNA[sec, ]
      #subb$C_Log2FC <- fold[core, ]$avg_log2FC
      #subb$S_Log2FC <- fold[sec, ]$avg_log2FC
      
      gl_res <- rbindlist(list(gl_res, subb))
    }
    
    fwrite(gl_res, paste0(dirr, "/final/SCOPE_Gene_Level.csv"))
    
    #GENE LEVEL

    genelevel <- fread(paste0(dirr, "/final/SCOPE_Gene_Level.csv"))
    genelevel[is.na(genelevel)] <- 0
    all_mgenes <- c()  
    results <- c()
    imp_pathways <- c()
    final <- c()

    
    for(x in nclss){
      print(c("Currently calculating for cluster ", x))
      sub_genelevel <- genelevel[genelevel$Cluster == x, ]

      #sub_genelevel <- sub_genelevel[sub_genelevel$`Is Differential Correlation Sig` == "Yes" | (is.na(sub_genelevel$`Other Correlation`) & sub_genelevel$`Cluster Correlation` >= diff_sig) | (is.na(sub_genelevel$`Other Correlation`) & sub_genelevel$`Cluster Correlation` <= -diff_sig) | (is.na(sub_genelevel$`Cluster Correlation`) & sub_genelevel$`Other Correlation` >=diff_sig)| (is.na(sub_genelevel$`Cluster Correlation`) & sub_genelevel$`Other Correlation` <= -diff_sig), ]
      
      
      if(nrow(sub_genelevel) == 0) {
        next}
      
      sub_genelevel$diffmet <- mapply(max,sub_genelevel$`Correlation`, sub_genelevel$`Cluster Correlation`, sub_genelevel$`Other Correlation`, sub_genelevel$`Differential Correlation`)
      
      pathway <- subset(sub_genelevel, select = c("Cluster", "Core Gene", "Secondary Gene", "KEGG GeneSet"))
      pathway_result <- data.frame(table(pathway$`KEGG GeneSet`))
      colnames(pathway_result) <- c("KEGG GeneSet", "Frequency")
      pathway_result$cluster <- x
      pathway_result <- merge(pathway_result, distinct(sub_genelevel[, c("KEGG GeneSet", "Pathway Name")]), by = "KEGG GeneSet", all = TRUE)
      pathway_result <- pathway_result[order(pathway_result$Freq, decreasing = TRUE), ]
      print("Important Pathways done")
      imp_pathways <- rbindlist(list(imp_pathways, pathway_result), fill = TRUE)
      
      sub_genelevel <- subset(sub_genelevel, select = c("Core Gene", "Secondary Gene", "KEGG GeneSet", "diffmet")) 
      diff_coregenes <- unique(sub_genelevel$`Core Gene`)
      final_results <- c()
      
      for(k in 1: length(diff_coregenes)){
        tmp2 <- sub_genelevel[sub_genelevel$`Core Gene` == diff_coregenes[[k]], ]
        ###Find if one core-secondary gene pair is involved in many pathways. 
        tmp3 <- subset(tmp2, select = c(`Core Gene`, `Secondary Gene`))
        tmp3 <- aggregate(list(numdup=rep(1,nrow(tmp3))), tmp3, length)
        tmp3 <- tmp3[order(tmp3$numdup, decreasing = TRUE), ]
        tmp2 <- distinct(subset(tmp2, select = -`KEGG GeneSet`))
        
        curr_result <- merge(tmp2, tmp3, by = c("Core Gene", "Secondary Gene"), all = TRUE)
        curr_result <- curr_result[order(curr_result$`diffmet`, decreasing= TRUE),]
        final_results <- rbindlist(list(final_results, curr_result), fill = TRUE)}
      
      final_results$cluster <- x
      
      
      npathways <- final_results$numdup
      npathways <- na.omit(npathways)
      max_pathways <- max(npathways)
      
      diff <- final_results$diffmet
      max_diff <- max(diff)
      
      final_results$adjusted_differential <- (final_results$diffmet / max_diff )* (final_results$numdup / max_pathways)
      adiff <- final_results$adjusted_differential
      adiff <- na.omit(adiff)
      max_adiff <- max(adiff)
      adj_cutoff <- 0.2 * max_adiff
      final_results <- final_results[final_results$adjusted_differential >= adj_cutoff, ]
      
      print("Adjusted Differential done")
      #final_results <- final_results[(final_results$adjusted_differential >= quantile(final_results$adjusted_differential,prob=1-adj_cutoff/100, na.rm = TRUE) | final_results$`Differential Correlation` >= quantile(final_results$`Differential Correlation`,prob=1-dif_cutoff/100, na.rm=TRUE)), ]
      final <- rbindlist(list(final_results, final), fill = TRUE)
      marker_genes = c()
      row_no = 1
      while(TRUE){
        genes <- c(final_results[row_no, ]$`Core Gene`, final_results[row_no, ]$`Secondary Gene`)
        marker_genes <- unique(c(marker_genes, genes))
        if(row_no ==nrow(final_results)) {
          break
        } else {
          row_no = row_no +1
        }
      }
      
      print("Calculating Fold Change")
      all_marker <- FoldChange(all_cell, ident.1 = x, features = unique(marker_genes)) ##############################################################################
      all_marker$Gene <- rownames(all_marker)
      all_marker$cluster = x
      all_mgenes <- rbindlist(list(all_mgenes, all_marker), fill = TRUE)
      marker <- all_marker[all_marker$pct.1 > 0.45 | all_marker$pct.2 > 0.45, ]
      marker <- marker[marker$avg_log2FC >= FC_thres | marker$avg_log2FC <= -FC_thres, ]
      results <- rbindlist(list(marker, results), fill = TRUE)
    }
    
    fwrite(all_mgenes, paste0(dirr, "/final/All_Marker_Genes.csv"))
    fwrite(results, paste0(dirr, "/final/Filtered_Marker_Genes.csv"))
    fwrite(imp_pathways, paste0(dirr, "/final/Important_Pathways.csv"))
    
    
    
    #results <- fread(paste0(dirr, "./final/Filtered_Marker_Genes.csv"))
    m_genes <- results$Gene
    f_data <- genelevel[genelevel$`Core Gene` %in% m_genes | genelevel$`Secondary Gene` %in% m_genes, ]
    
    f_data <- subset(f_data, select = c(Cluster, `Core Gene`, `Secondary Gene`, `KEGG GeneSet`, `Pathway Name`))
    mel_data <- melt(f_data, id.vars = c("Cluster", "KEGG GeneSet", "Pathway Name"))
    
    ncls = unique(mel_data$Cluster)
    
    fin_res = c()
    for(i in ncls){
      
      m_genes2 <- results[results$cluster == i,]$Gene
      mel_data2 <- mel_data[mel_data$Cluster == i & mel_data$value %in% m_genes2,]
      f_table2 <- data.frame(table(mel_data2$value, mel_data2$`Pathway Name`))
      un_genes <- unique(f_table2$Var1)
      for(gene in un_genes){
        
        m1 <- f_table2[f_table2$Var1 == gene, ]
        m1 <- m1[m1$Freq !=0, ]
        t <- paste(as.character(m1[order(m1$Freq, decreasing = TRUE),][1:10,]$Var2), collapse = ",")
        t <- gsub(",NA", "", t)
        cur_res = list("Gene" = gene, "pathways" = t, "cluster" = i)
        fin_res = rbindlist(list(fin_res, cur_res))
      }
    }
    fwrite(fin_res, paste0(dirr, "/final/path_gene.csv"))
  }
}
