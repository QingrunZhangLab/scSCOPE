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

diffCorrelation_null <- function(ex.data, cls1, cls2 = NULL, START_SEED = 2222, PROBES_PER_ITER = 1000, ITERS = 100){
  
  if(!is.null(cls2)){
    expr.data <- ex.data[ex.data$phenotype %in% c(cls1, cls2), ]
    if(file.exists(paste0("diffCox/scope", cls1, "_", cls2, "_DiffCorrNull.csv"))){
      return()
    }
  }else{
    expr.data <- ex.data
    if(file.exists(paste0("diffCox/scope", cls1, "_DiffCorrNull.csv"))){
      return()
    }
  }
  expr.data$phenotype <- ifelse(expr.data$phenotype == cls1, 0,1)
  print("Current phenotypes after substituting: ")
  print(table(expr.data$phenotype))
  
  cluster1 <- expr.data[expr.data$phenotype == 0,]
  cluster2 <- expr.data[expr.data$phenotype == 1,]
  null_dist <- c()
  set.seed(START_SEED)
  for(i in 1:ITERS){
    print(paste0("Current trial ", i))
    rnd_smp <- sample(colnames(expr.data)[2:(ncol(expr.data)-1)], PROBES_PER_ITER)
  
    cluster1_sel <- cluster1[, ..rnd_smp]
    cluster2_sel <- cluster2[, ..rnd_smp]
    
    cluster1_cor <- cor(cluster1_sel)
    cluster2_cor <- cor(cluster2_sel)
    
    cluster1_corflat <- flattenCorrMatrix(cluster1_cor)
    cluster2_corflat <- flattenCorrMatrix(cluster2_cor)
    
    merged <- merge(cluster1_corflat, cluster2_corflat, by = c("row", "column"))
    merged$diffcor <- abs(merged$cor.x - merged$cor.y)
    
    if(is.null(null_dist)){
      null_dist <- as.data.table(merged$diffcor)
    }else{
      null_dist <- cbind(null_dist, merged$diffcor)
    }
  }
  print("Writing all results.")
  
  if(!is.null(cls2)){
    fwrite(null_dist, paste0("./diffCox/", PREFIX, cls1, "_", cls2, "_DiffCorrNull.csv"))
  }else{
    fwrite(null_dist, paste0("./diffCox/", PREFIX, cls1, "_DiffCorrNull.csv"))
  }
}

sparse_lasso <- function(ex.data, cls1, cls2 = NULL, ITERS_L = 200, FOLDS = 10){
  if(!is.null(cls2)){
    expr.data <- ex.data[ex.data$phenotype %in% c(cls1, cls2), ]
    if(file.exists(paste0("results_sparse/", PREFIX, cls1, "_", cls2, "_Summary.csv"))){
      return()
    }
  }else{
    if(file.exists(paste0("results_sparse/", PREFIX, cls1, "_Summary.csv"))){
      return()
    }
    expr.data <- ex.data
  }
  expr.data$phenotype <- ifelse(expr.data$phenotype == cls1, 0,1)
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
  
  
  if(!is.null(cls2)){
    fwrite(cof_sum, paste0("results_sparse/", PREFIX, cls1,"_", cls2,  "_Summary.csv"))
    fwrite(results, paste0("results_sparse/", PREFIX, cls1, "_", cls2, "_AllCoefs.csv"))
    fwrite(pred_metrics, paste0("results_sparse/", PREFIX,cls1, "_", cls2, "_PredMetrics.csv"))
  } else{
    fwrite(cof_sum, paste0("results_sparse/", PREFIX, cls1,  "_Summary.csv"))
    fwrite(results, paste0("results_sparse/", PREFIX, cls1, "_AllCoefs.csv"))
    fwrite(pred_metrics, paste0("results_sparse/", PREFIX, cls1, "_PredMetrics.csv"))
   
  }
}

coexpression <- function(ex.data, CORE_CUTOFF = 160, sub_sample = 0.6, corr_iter_cutoff = 80, iter_cor = 100,
                        CORR_PERCENTILE_THRESHOLD = 0.975, DIFFCORR_PERCENTILE_THRESHOLD = 0.975, cls1, cls2= NULL){
  
  expr.data <- ex.data
  if(!is.null(cls2)){
    SUMMARY_FILE = paste0(PREFIX, cls1, "_", cls2, "_Summary.csv")
    y = str_remove(SUMMARY_FILE, "_Summary.csv")
    expr.data <- expr.data[expr.data$phenotype %in% c(cls1, cls2), ]
    if(file.exists(paste0("results_sparse/corSparse/", y, "_SecondaryGenes.csv"))){
      return()
    }
  } else{
    SUMMARY_FILE = paste0(PREFIX, cls1, "_Summary.csv")
    y = str_remove(SUMMARY_FILE, "_Summary.csv")
    if(file.exists(paste0("results_sparse/corSparse/", y, "_SecondaryGenes.csv"))){
      return()
    }
  }
  
  CORR_NULL_FILE = paste0("correlation/", PREFIX, "_CorrNull.csv")
  
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
  colnames(gene_coords) <- "gene_name"
  
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
  
  print("Reading null distribution of differential co-expressions...")
  diffcorr_null_dist <- fread(DIFFCORR_NULL_FILE)
  diffcorr_cutoff <- c()
  print("Calculating median percentile from differential co-expression data...")
  for(i in 1:ncol(diffcorr_null_dist)){
    diffcorr_cutoff <- c(diffcorr_cutoff, quantile(diffcorr_null_dist[, ..i][[1]], DIFFCORR_PERCENTILE_THRESHOLD, na.rm = TRUE)[[1]])
  }
  
  diffcorr_cut <- median(diffcorr_cutoff)

  
  cluster0 <- expr.data[expr.data$phenotype == cls1,]
  other0 <- expr.data[expr.data$phenotype !=cls1,]
  
  corr_res = c()
  
  for(iter in 1:iter_cor){
    
    cluster <- stratified(cluster0, "phenotype", sub_sample )
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
      corr_res <- merge(corr_res, curr_final, by = c("Core Gene", "Secondary Gene"), all.x=TRUE, all.y=TRUE,)
      
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
  fwrite(results, paste0("./results_sparse/corSparse/", y, "_AllCoreCorrelations.csv"))
  fwrite(filtered_results, paste0("./results_sparse/corSparse/", y, "_SecondaryGenes.csv"))
  
}


pathwayEnrichment <- function(pathwayDatabase = c("KEGG", "GOBP", "GOCC", "GOMF"), pathgeneId, organism, cls1, cls2 = NULL){
  pdb <- unique(pathwayDatabase)
  print(cls2)
  if(!is.null(cls2)){
    SECONDARY_GENE_FILE = paste0(PREFIX, cls1, "_", cls2, "_SecondaryGenes.csv")
  } else{
    SECONDARY_GENE_FILE = paste0(PREFIX, cls1, "_SecondaryGenes.csv")
  }
  x = SECONDARY_GENE_FILE
  full_data <- fread(paste0("results_sparse/corSparse/", SECONDARY_GENE_FILE))
  
  print("File Read")
  
  if("KEGG" %in% pdb){
    if(exists("final_path")){
      rm(final_path)}
    if(!file.exists(paste0("results_sparse/corSparse/all_path/path/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))){
      dir.create("results_sparse/corSparse/all_path/path")
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
          path_res$cluster <- str_remove(x, "_Summary.csv")
          
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
      dir.create("results_sparse/corSparse/all_path/cc")
      
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
          path_res$cluster <- str_remove(x, "_Summary.csv")
          
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
      dir.create("results_sparse/corSparse/all_path/mf")
      
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
          path_res$cluster <- str_remove(x, "_Summary.csv")
          
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
      dir.create("results_sparse/corSparse/all_path/bp")
      
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
          path_res$cluster <- str_remove(x, "_Summary.csv")
          
          if(exists("final_bp")){
            final_bp <- rbindlist(list(final_bp, path_res))
          } else{
            final_bp <- path_res
          }
        }
      }
    }
    
  }
  print(x)
  if(exists("final_bp")){
    fwrite(final_bp, paste0("results_sparse/corSparse/all_path/bp/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))}
  
  if(exists("final_cc")){
    fwrite(final_cc, paste0("results_sparse/corSparse/all_path/cc/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))}
  
  if(exists("final_mf")){
    fwrite(final_mf, paste0("results_sparse/corSparse/all_path/mf/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))}
  
  if(exists("final_path")){
    fwrite(final_path, paste0("results_sparse/corSparse/all_path/path/", str_remove(x, "_SecondaryGenes.csv"), "_Pathways.csv"))}
}



markergenes <- function(ex.data, DIFFCORR_PERCENTILE_THRESHOLD = 0.975, FC_thres = 0.5, cls1, cls2 = NULL){
  expr.data <- data.frame(as.matrix(ex.data, rownames = 1))
  for(dirr in list.dirs(path = "./results_sparse/corSparse/all_path", recursive = FALSE)){
    dir.create(paste0(dirr, "/final"), showWarnings = FALSE)
    if(!is.null(cls2)){
      all_pathways <- fread(paste0(dirr, "/", PREFIX, cls1, "_", cls2, "_Pathways.csv"))
    }else{
      all_pathways <- fread(paste0(dirr, "/", PREFIX, cls1,"_Pathways.csv"))
    }
    all_pathways$cluster <- str_replace(all_pathways$cluster, "_SecondaryGenes.csv", "")
    all_pathways$cluster <- str_replace(all_pathways$cluster, PREFIX, "")
      
    all_corrs = c()
    
    if(!is.null(cls2)){
      cls_no <- paste0(cls1, "_", cls2)
      all_corrs <- fread(paste0( "./results_sparse/corSparse/", PREFIX, cls1, "_", cls2, "_SecondaryGenes.csv"))
      diffcorr_null_dist <- fread(paste0("./diffCox/", PREFIX, cls1, "_", cls2, "_DiffCorrNull.csv"))
    }else{
      cls_no = cls1
      all_corrs <- fread(paste0( "./results_sparse/corSparse/", PREFIX, cls1, "_SecondaryGenes.csv"))
      diffcorr_null_dist <- fread(paste0("./diffCox/", PREFIX, cls1, "_DiffCorrNull.csv"))
    }
    
    diffcorr_cutoff <- c()
    print("Calculating median percentile from differential co-expression data...")
    for(i in 1:ncol(diffcorr_null_dist)){
      diffcorr_cutoff <- c(diffcorr_cutoff, quantile(diffcorr_null_dist[, ..i][[1]], DIFFCORR_PERCENTILE_THRESHOLD, na.rm = TRUE)[[1]])
    }
    diffcorr_cut <- median(diffcorr_cutoff)
    all_corrs$`Differential Correlation Sig` <- diffcorr_cut
    all_corrs$`Is Differential Correlation Sig` <- ifelse(all_corrs$`Differential Correlation` > all_corrs$`Differential Correlation Sig`, "Yes", "No")
    all_corrs$cluster = as.character(cls_no)
      
    final_data <- c()
    for(i in 1:nrow(all_pathways)){
      genes <- unique(str_split(all_pathways$userId[i], ";")[[1]])
      curr_data <- data.table(cluster = rep(all_pathways$cluster[i], length(genes)),
                              `KEGG GeneSet` = rep(all_pathways$geneSet[i], length(genes)),
                              `Pathway Name` = rep(all_pathways$description[i], length(genes)),
                              `Core Gene` = rep(all_pathways$`Core Gene`[i], length(genes)),
                              `Secondary Gene` = genes)
      
      final_data <- rbindlist(list(final_data, curr_data))
    }
    
    final_data$`Core Gene` <- all_corrs$`Core Gene`[match(final_data$`Core Gene`, all_corrs$`Core Gene`)]
    final_data$`Secondary Gene` <- all_corrs$`Secondary Gene`[match(final_data$`Secondary Gene`, all_corrs$`Secondary Gene`)]
    final_data$cluster <- str_remove(final_data$cluster, "_SecondaryGenes.csv")
    
    final_data <- merge(final_data, all_corrs, by = c("cluster", "Core Gene", "Secondary Gene"), no.dups = TRUE)
    fwrite(final_data, paste0(dirr, "/final/SCOPE_Gene_Level", cls_no, ".csv"))
    
    #GENE LEVEL
    genelevel <- fread(paste0(dirr, "/final/SCOPE_Gene_Level", cls_no, ".csv"))
    genelevel[is.na(genelevel)] <- 0
    
    all_mgenes <- c()  
    results <- c()
    final <- c()

    genelevel$diffmet <- mapply(max, genelevel$`Correlation`, genelevel$`Cluster Correlation`, genelevel$`Other Correlation`, genelevel$`Differential Correlation`)
    
    pathway <- subset(genelevel, select = c("cluster", "Core Gene", "Secondary Gene", "KEGG GeneSet"))
    imp_pathways <- data.frame(table(pathway$`KEGG GeneSet`))
    colnames(imp_pathways) <- c("geneSet", "Frequency")
    imp_pathways$cluster <- cls_no
    imp_pathways <- merge(imp_pathways, distinct(all_pathways[, c("geneSet", "description")]), by = "geneSet", all = TRUE)
    imp_pathways <- imp_pathways[order(imp_pathways$Freq, decreasing = TRUE), ]
    print("Important Pathways done")

    genelevel <- subset(genelevel, select = c("Core Gene", "Secondary Gene", "KEGG GeneSet", "diffmet", "Pathway Name", "cluster")) 
    diff_coregenes <- unique(genelevel$`Core Gene`)
    final_results <- c()
    
    for(k in 1: length(diff_coregenes)){
      tmp2 <- genelevel[genelevel$`Core Gene` == diff_coregenes[[k]], ]
      ###Find if one core-secondary gene pair is involved in many pathways. 
      tmp3 <- subset(tmp2, select = c(`Core Gene`, `Secondary Gene`))
      tmp3 <- aggregate(list(numdup=rep(1,nrow(tmp3))), tmp3, length)
      tmp3 <- tmp3[order(tmp3$numdup, decreasing = TRUE), ]
      tmp2 <- distinct(subset(tmp2, select = -`KEGG GeneSet`))
      
      curr_result <- merge(tmp2, tmp3, by = c("Core Gene", "Secondary Gene"), all = TRUE)
      curr_result <- curr_result[order(curr_result$`diffmet`, decreasing= TRUE),]
      final_results <- rbindlist(list(final_results, curr_result), fill = TRUE)}
    
    final_results$cluster <- cls_no
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
    
    print("Fold change doing")
    ident1 <- str_split(cls_no, "_")[[1]][1]
    ident2 <- str_split(cls_no, "_")[[1]][2]
    if(is.na(ident2)){
      ident2 = NULL
    }
    all_marker <- fold_change(ex.data, cls1 = ident1, cls2 = ident2, features = unique(marker_genes)) ##############################################################################
    all_marker$Gene <- rownames(all_marker)
    all_marker$cluster = cls_no
    all_mgenes <- rbindlist(list(all_mgenes, all_marker), fill = TRUE)
    marker <- all_marker[all_marker$pct.1 > 0.45 | all_marker$pct.2 > 0.45, ]
    marker <- marker[marker$Avg_log2FC >= FC_thres | marker$Avg_log2FC <= -FC_thres, ]
    
    #mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    #gene_coords <- data.frames$GeneName <- gene_coords$external_gene_name
  }
  
  fwrite(all_mgenes, paste0(dirr, "/final/All_Marker_Genes", cls_no, ".csv"))
  fwrite(marker, paste0(dirr, "/final/Filtered_Marker_Genes", cls_no, ".csv"))
  fwrite(imp_pathways, paste0(dirr, "/final/Important_Pathways.csv", cls_no, ".csv"))
  
  
  
  #results <- fread(paste0(dirr, "./final/Filtered_Marker_Genes.csv"))
  m_genes <- marker$Gene
  f_data <- genelevel[genelevel$`Core Gene` %in% m_genes | genelevel$`Secondary Gene` %in% m_genes, ]
  
  
  
  f_data <- subset(f_data, select = c(cluster, `Core Gene`, `Secondary Gene`, `KEGG GeneSet`, `Pathway Name`))
  mel_data <- melt(f_data, id.vars = c("cluster", "KEGG GeneSet", "Pathway Name"))
  mel_data$cluster <- str_replace(mel_data$cluster, PREFIX, "")

  m_genes2 <- marker[marker$cluster == cls_no,]$Gene
  mel_data2 <- mel_data[mel_data$cluster == cls_no & mel_data$value %in% m_genes2,]
  f_table2 <- data.frame(table(mel_data2$value, mel_data2$`Pathway Name`))
  un_genes <- unique(f_table2$Var1)
  fin_res = c()
  for(gene in un_genes){
    m1 <- f_table2[f_table2$Var1 == gene, ]
    m1 <- m1[m1$Freq !=0, ]
    t <- paste(as.character(m1[order(m1$Freq, decreasing = TRUE),][1:10,]$Var2), collapse = ",")
    t <- gsub(",NA", "", t)
    cur_res = list("Gene" = gene, "pathways" = t, "cluster" = cls_no)
    fin_res = rbindlist(list(fin_res, cur_res))
  }
  fwrite(fin_res, paste0(dirr, "/final/path_gene", cls_no, ".csv"))
}




fold_change <- function(ex.data, cls1, cls2 = NULL, features){
  thresh.min <- 0
  cells.1 <- ex.data[ex.data$phenotype == cls1, ]$V1
  
  if(!(is.null(cls2))){
    cells.2 <- ex.data[ex.data$phenotype == cls2, ]$V1
  }else{
    cells.2 <- ex.data[ex.data$phenotype != cls1, ]$V1
  }
  ex.data2 <- data.frame(ex.data)
  rownames(ex.data2) <- ex.data$V1
  ex.data2 <- ex.data2[, -1]
  ex.data2 <- t(ex.data2)
  thresh.min <- 0
  grp1 <- ex.data2[features, cells.1, drop = FALSE]
  grp2 <- ex.data2[features, cells.2, drop = FALSE]
  
  pct.1 <- round(
    x = rowSums(x = grp1 > thresh.min) /
      length(x = cells.1),
    digits = 3
  )
  pct.2 <- round(
    x = rowSums(x = grp2 > thresh.min) /
      length(x = cells.2),
    digits = 3
  )
  log2FC <- rowMeans((grp1) - rowMeans(grp2))
  
  results <- data.frame("Avg_log2FC" = log2FC, "pct.1" = pct.1, "pct.2" = pct.2)
  return(results)
} 




findMarkerGenes <- function(ex.data, cluster1, cluster2=NULL, corr_percentile_threshold = 0.975, 
                            probes_per_iter = 1000, iters = 100, iters_lasso = 200, folds = 10, 
                            start_seed = 2222, org, pathDB, core_cutoff = 160, iter_correlation = 100, 
                            corr_iter_cut = 80, geneId,diff_threshold = 0.975, FC_threshold = 0.5 )
{
  correlation_null(ex.data, ITERS = iters, CORR_PERCENTILE_THRESHOLD = corr_percentile_threshold, 
                   PROBES_PER_ITER= probes_per_iter)
  diffCorrelation_null(ex.data, cls1 = cluster1, cls2=cluster2,START_SEED = start_seed,
                       PROBES_PER_ITER= probes_per_iter, ITERS = iters)
  sparse_lasso(ex.data, cls1 = cluster1, cls2 = cluster2,ITERS_L = iters_lasso, FOLDS = folds)
  
  coexpression(ex.data, CORE_CUTOFF = core_cutoff, iter_cor = iter_correlation, 
               corr_iter_cutoff = corr_iter_cut, cls1 = cluster1, cls2 = cluster2)
  
  pathwayEnrichment(pathgeneId = geneId, organism = org, pathwayDatabase = pathDB, cls1 = cluster1,
                    cls2 = cluster2)
  markergenes(ex.data, DIFFCORR_PERCENTILE_THRESHOLD = diff_threshold, cls1 = cluster1, 
              FC_thres = FC_threshold, cls2 = cluster2)
}


findAllMarkerGenes <- function(ex.data, cluster1, cluster2=NULL, corr_percentile_threshold = 0.975, 
                               probes_per_iter = 1000, iters = 100, iters_lasso = 200, folds = 10, 
                               start_seed = 2222, org, pathDB, core_cutoff = 160, iter_correlation = 100, 
                               corr_iter_cut = 80, geneId,diff_threshold = 0.975, FC_threshold = 0.5 )
{
  for(i in unique(ex.data$phenotype)){
    correlation_null(ex.data, ITERS = iters, CORR_PERCENTILE_THRESHOLD = corr_percentile_threshold, 
                     PROBES_PER_ITER= probes_per_iter)
    diffCorrelation_null(ex.data, cls1 = cluster1, cls2=cluster2,START_SEED = start_seed,
                         PROBES_PER_ITER= probes_per_iter, ITERS = iters)
    sparse_lasso(ex.data, cls1 = cluster1, cls2 = cluster2,ITERS_L = iters_lasso, FOLDS = folds)
    
    coexpression(ex.data, CORE_CUTOFF = core_cutoff, iter_cor = iter_correlation, 
                 corr_iter_cutoff = corr_iter_cut, cls1 = cluster1, cls2 = cluster2)
    
    pathwayEnrichment(pathgeneId = geneId, organism = org, pathwayDatabase = pathDB, cls1 = cluster1,
                      cls2 = cluster2)
    markergenes(ex.data, DIFFCORR_PERCENTILE_THRESHOLD = diff_threshold, cls1 = cluster1, 
                FC_thres = FC_threshold, cls2 = cluster2)
  }
  
  
}

findAllMarkerGenes1vs1 <- function(ex.data, cluster1, cluster2=NULL, corr_percentile_threshold = 0.975, 
                                   probes_per_iter = 1000, iters = 100, iters_lasso = 200, folds = 10, 
                                   start_seed = 2222, org, pathDB, core_cutoff = 160, iter_correlation = 100, 
                                   corr_iter_cut = 80, geneId,diff_threshold = 0.975, FC_threshold = 0.5 ){
  doneClusters <- c()
  for(cls1 in unique(ex.data$phenotype)){
    for(cls2 in unique(ex.data$phenotype)){
      if(cls1 != cls2 & !(cls2 %in% doneClusters)){
        correlation_null(ex.data, ITERS = iters, CORR_PERCENTILE_THRESHOLD = corr_percentile_threshold, 
                         PROBES_PER_ITER= probes_per_iter)
        diffCorrelation_null(ex.data, cls1 = cluster1, cls2=cluster2,START_SEED = start_seed,
                             PROBES_PER_ITER= probes_per_iter, ITERS = iters)
        sparse_lasso(ex.data, cls1 = cluster1, cls2 = cluster2,ITERS_L = iters_lasso, FOLDS = folds)
        
        coexpression(ex.data, CORE_CUTOFF = core_cutoff, iter_cor = iter_correlation, 
                     corr_iter_cutoff = corr_iter_cut, cls1 = cluster1, cls2 = cluster2)
        
        pathwayEnrichment(pathgeneId = geneId, organism = org, pathwayDatabase = pathDB, cls1 = cluster1,
                          cls2 = cluster2)
        markergenes(ex.data, DIFFCORR_PERCENTILE_THRESHOLD = diff_threshold, cls1 = cluster1, 
                    FC_thres = FC_threshold, cls2 = cluster2)
      }
    }
    doneClusters <- c(doneClusters, cls1)
    
  }
  
}


mgenes_correlation <- function(ex.data, cls1, cls2 = NULL, mgenes ){
    if(!is.null(cls2)){
      expr.data <- ex.data[ex.data$phenotype %in% c(cls1, cls2), ]
    }else{
      expr.data <- ex.data
    }
		expr.data$phenotype <- ifelse(expr.data$phenotype == cls1, 0,1)
		cluster <- expr.data[expr.data$phenotype == 0,]
		other <- expr.data[expr.data$phenotype == 1,]
		results <- c()
		
		clsM <- as.matrix(cluster[, !c("V1", "phenotype")])
		otherM <- as.matrix(other[, !c("V1", "phenotype")])
		allM <- as.matrix(expr.data[,!c("V1", "phenotype")])
		
		clsM_core <- as.matrix(subset(cluster, select = mgenes))
		otherM_core <- as.matrix(subset(other, select = mgenes))
		allM_core <- as.matrix(subset(expr.data, select = mgenes))
		
		
		cls_res <- corSparse(clsM_core, clsM)
		other_res <- corSparse(otherM_core, otherM)
		
		rownames(cls_res) <- colnames(clsM_core)
		colnames(cls_res) <- colnames(clsM)
		
		rownames(other_res) <- colnames(otherM_core)
		colnames(other_res) <- colnames(otherM)
		
		results_cluster <- data.table(`Core Gene` = rownames(cls_res), cls_res)
		results_cluster <- melt(results_cluster, id.vars = "Core Gene", value.name = "Cluster Correlation", variable.name = "Secondary Gene")
		
		results_other <- data.table(`Core Gene` = rownames(other_res), other_res)
		results_other <- melt(results_other, id.vars = "Core Gene", value.name = "Other Correlation", variable.name = "Secondary Gene")

		results <- merge(results_cluster, results_other, by = c("Core Gene", "Secondary Gene"), all = TRUE)	
	
	  results[is.na(results)] <- 0
  	return(results)
}
plotData <- function(all_corr, s_path){
  all_corr$`Core Gene Name` <- all_corr$`Core Gene`
  all_corr$`Secondary Gene Name` <- all_corr$`Secondary Gene`
  scope_plot = c()
  for(curr_pathway in s_path){
    genes = c()
    count = 0
    while(length(genes)==0 && count < 50){
      query <- tryCatch(keggGet(curr_pathway), error = function(e){NULL})
      genes <- query[[1]]$GENE
      print(c(length(genes), count))
      count <- count+1
    }
    
    if(is.null(genes)){
      print(c("No genes identified for the pathway", curr_pathway))
      next
    }
    print(c("KEGG Query done", curr_pathway))
    genes <- genes[seq(2, length(genes), 2)]
    genes <- str_split(genes, ";", simplify = TRUE)[,1]
    
    sub_corr <- all_corr[toupper(all_corr$`Secondary Gene Name`) %in% toupper(genes), ]
    sub_corr$`KEGG GeneSet` <- curr_pathway 
    
    col_order <- c("KEGG GeneSet", "Core Gene", "Core Gene Name", "Secondary Gene", "Secondary Gene Name", "Cluster Correlation", "Other Correlation")
    sub_corr <- subset(sub_corr, select = col_order)
    
    scope_plot <- rbindlist(list(sub_corr, scope_plot), fill = TRUE)}
  return(scope_plot)
}
extend_with_bootstrap <- function(data1, data2, desired_length) {
  set.seed(123)
  original_length <- length(data1)
  additional_samples_needed <- desired_length - original_length
  if (additional_samples_needed > 0) {
    additional_indices <- sample(1:original_length, additional_samples_needed, replace = TRUE)
    extended_data1 <- c(data1, data1[additional_indices])
    extended_data2 <- c(data2, data2[additional_indices])
    
  } else {
    extended_data1 <- data1
    extended_data2 <- data2
  }
  return(list(extended_data1, extended_data2))
}
plot_file <- function(plot_data, mgene){
  curr_index <- 1
  plot_collection <- list()
  pathh <- unique(plot_data$`KEGG GeneSet`)
  print(length(pathh))
  diff_path = c()
  for(path in pathh){
    curr_plot_data <- plot_data[`KEGG GeneSet` == path]
    tests <- c()
    tmpx <- curr_plot_data$`Cluster Correlation`
    tmpy <- curr_plot_data$`Other Correlation`
    
    ress <- extend_with_bootstrap(tmpx, tmpy, 100)
    tmpx <- ress[[1]]
    tmpy <- ress[[2]]
    tests <- data.table(cluster = c,
                        group1 = "Cluster",
                        group2 = "Other",
                        p = signif(ks.test(tmpx, tmpy)$p.value))
    
    diff_path = append(diff_path, path)
    sec_genes <- unique(curr_plot_data$`Secondary Gene Name`)
    
    for(t in c("Cluster", "Other")){
      curr_plot_data <- curr_plot_data[`KEGG GeneSet` == path,]
      curr_nodes_data <- unique(curr_plot_data[, c("Core Gene Name", "Secondary Gene Name")])
      curr_nodes_data <- melt(curr_nodes_data, variable.name = "type", value.name = "name", measure.vars = c("Core Gene Name", 
                                                                                                             "Secondary Gene Name"))
      curr_nodes_data$type <- str_split(curr_nodes_data$type, pattern = " ", simplify = TRUE)[,1]
      curr_nodes_data <- unique(curr_nodes_data)
      curr_nodes_data <- setcolorder(curr_nodes_data, c("name", "type"))
      curr_nodes_data$scale <- ifelse(curr_nodes_data$type == "Core", 10, 2)
      
      if(t == "Cluster"){
        curr_edges <- curr_plot_data[, c("Core Gene Name", "Secondary Gene Name", "Cluster Correlation")]
      }else if(t == "Other")
      {
        curr_edges <- curr_plot_data[, c("Core Gene Name", "Secondary Gene Name", "Other Correlation")]
      }
      
      colnames(curr_edges) <- c("from", "to", "corr")
      curr_graph <- tbl_graph(nodes = curr_nodes_data, edges = curr_edges)
      p <- ggraph(curr_graph, layout = "star") +
        geom_edge_link(aes(colour = corr), edge_width = 0.5, show.legend = c(colour = TRUE, edge_width = FALSE))+
        # scale_edge_colour_gradientn(colours = c("red", "white","blue"),
        #                             breaks = c(-1,0,1), limits=c(-1,1),
        #                             guide = "edge_colourbar")+
        scale_edge_colour_gradient2(
          limits = c(-1,1),
          low = "red",
          mid = "white",
          high = "blue",
          guide = "edge_colourbar",
        ) +
        xlim(-1.2, 1.2) + 
        ylim(-1.2, 1.2) +
        geom_node_point(aes(colour = type, size = scale), show.legend = FALSE) +
        scale_color_manual(values = c("Core" = "skyblue", "Secondary" = "grey")) +
        geom_node_text(aes(label=name, 
                           filter = type=="Core", 
                           fontface = "bold"),
                       size = 9/.pt,
                       show.legend = FALSE) +
        # geom_node_text(aes(label=name,
        #                    filter = type=="Secondary",
        #                    fontface = "bold"),
        #                size = 4/.pt,
        #                show.legend = FALSE) +
        coord_cartesian(clip = "off") +
        theme_void() +
        theme(text = element_text(size = 9),
              legend.title = element_blank(),
              legend.text = element_text(size = 6),
              legend.title.align = 0.5,
              legend.key.height= unit(0.3, 'cm'),
              legend.key.width= unit(0.3, 'cm'),
              plot.margin = margin(t = 0, r = -1, b = 0, l = -1, unit = "cm"),
              panel.spacing = unit(-1, units = "cm"))
      
      if(curr_index == 1){
        p <- p + annotation_custom(grob = textGrob(label = "Cluster Cor.", 
                                                   hjust = "center",
                                                   rot = 90,
                                                   gp = gpar(cex = .6)),
                                   xmax = -2)
        
        p <- p + labs(title = stringr::str_wrap(pathh, width=15)) +
          theme(plot.title = element_text(face = "bold", size = 7, hjust = 0.5))
        
        
      }
      
      if(curr_index == 2){
        p <- p + annotation_custom(grob = textGrob(label = "Other Cor.", 
                                                   hjust = "center",
                                                   rot = 90,
                                                   gp = gpar(cex = .6)),
                                   xmax = -2)
      }
      
      if(curr_index%%3 ==1){
        p <- p + labs(title = stringr::str_wrap(pathh, width=15)) +
          theme(plot.title = element_text(face = "bold", size = 7, hjust = 0.5))
      }	 
      
      plot_collection[[curr_index]] <- p
      curr_index <- curr_index + 1
    }
    
    boxplot_data <- melt(curr_plot_data[, ], id.vars = c("Secondary Gene Name"), measure.vars = c("Cluster Correlation", "Other Correlation"))
    colnames(boxplot_data)[2:3] <- c("Tissue", "Correlation")
    boxplot_data$Tissue <- str_remove(boxplot_data$Tissue, " Correlation")
    
    bxplot <- ggboxplot(boxplot_data, x = "Tissue", y = "Correlation", color = "Tissue",
                        xlab = FALSE, legend = "right", ggtheme = theme_light(),
                        legend.title = "",
                        legend.position = "none",
                        palette = "npg",
                        order = c("Cluster", "Other"),
                        size = 0.05,
                        width = 0.6)+
      stat_pvalue_manual(tests,
                         label.size = 2,
                         y.position = 1.2,
                         label = "p = {p}") +
      scale_y_continuous(limits = c(-1, 1.3)) +
      theme(axis.line = element_line(),
            text = element_text(size = 7),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "none")
    
    if(curr_index != 3){
      bxplot <- bxplot + theme(axis.text.y = element_blank(),
                               axis.title.y = element_blank())
    }
    
    plot_collection[[curr_index]] <- bxplot
    curr_index <- curr_index + 1
  }
  
  pdiff <- wrap_plots(plot_collection, nrow = 3, ncol = length(diff_path)+1, byrow = FALSE, guides = "collect",
                      widths = c(rep(20, length(diff_path)), 2))
  
  p <- pdiff 
  plot_width = 3*(length(diff_path)) + 3
  return(p)
}

geneNetworkPlots <- function(expr.data, marker_gene, pathways_of_interest, cluster1, cluster2){
  expr.data <- dataFiltering(expr.data)
  mgene_corr_matrix <- mgenes_correlation(expr.data, cluster1, cluster2, marker_gene)
  
  all_plots <- c()
  for(i in marker_gene){
    mgene_corr_matrix_sep <- mgene_corr_matrix[mgene_corr_matrix$`Core Gene` == i, ]
    plot_data <- plotData(mgene_corr_matrix_sep, pathways_of_interest)
    plots <- plot_file(plot_data, i)
    all_plots[[i]] <- plots
  }
  return(all_plots)
}


diffCorr <- function(cluster, other, en_gen){
  clsM <- as.matrix(subset(cluster, select = en_gen))
  otherM <- as.matrix(subset(other, select = en_gen))
  
  cls_res <- corSparse(clsM, clsM)
  oth_res <- corSparse(otherM, otherM)
  
  rownames(cls_res) <- colnames(clsM)
  colnames(cls_res) <- colnames(clsM)
  
  rownames(oth_res) <- colnames(otherM)
  colnames(oth_res) <- colnames(otherM)
  
  res_cls <- data.table("from" = rownames(cls_res), cls_res)
  res_cls <- melt(res_cls, id.vars = "from", value.name = "weights", variable.name = "to")
  
  res_cls[is.na(res_cls$weights)]$weights <- 0
  
  
  res_oth <- data.table("from" = rownames(oth_res), oth_res)
  res_oth <- melt(res_oth, id.vars = "from", value.name = "weights", variable.name = "to")
  
  
  res_oth[is.na(res_oth$weights)]$weights <- 0
  
  res_cls <- res_cls[res_cls$from != res_cls$to, ]
  res_oth <- res_oth[res_oth$from != res_oth$to, ]
  tmps <- res_cls$`weights` - res_oth$`weights`
  return(tmps)
}

diffExp <- function(cluster, other, en_gen){
  clsM <- as.matrix(subset(cluster, select = en_gen))
  otherM <- as.matrix(subset(other, select = en_gen))
  nodes <- data.frame(id = en_gen)
  nodes$cls_exp <- sapply(nodes$id, function(x) mean(clsM[,x]))
  nodes$oth_exp <- sapply(nodes$id, function(x) mean(otherM[,x]))
  nodes$diff <- nodes$cls_exp - nodes$oth_exp
  return(nodes$diff)
}

plots <- function(res, nodes, c, cluster){
  links <- res
  links <- links[!(links$from == links$to), ]
  links <- links[!is.na(links$`weights`), ]
  links <- links[order(links$weights), ]
  links <- links %>% filter(row_number() %%2 == 1)
  links <- links[abs(links$weights) >= 0.2,]
  vis.nodes <- nodes
  vis.links <- links
  
  vis.nodes$shape  <- "dot"  
  vis.nodes$shadow <- TRUE # Nodes will drop shadow
  vis.nodes$title  <- vis.nodes$id # Text on click
  vis.nodes$label  <- vis.nodes$id # Node label
  vis.nodes$size   <- vis.nodes$cls_exp * 30 # Node size
  vis.nodes$borderWidth <- 1 # Node border width
  
  vis.nodes$color.background <- c("slategrey", "tomato")[as.integer(nodes$type == "MARKER") + 1]
  vis.nodes$color.border <- "black"
  vis.nodes$color.highlight.background <- "orange"
  vis.nodes$color.highlight.border <- "darkred"
  
  
  vis.links$width <- links$weights * 10 # line width
  vis.links$color <- c("red", "blue")[as.integer(links$weights >= 0) + 1]    # line color   # arrows: 'from', 'to', or 'middle'
  vis.links$smooth <- FALSE    # should the edges be curved?
  vis.links$shadow <- FALSE    # edge shadow
  
  #vis.links$width <- abs(vis.links$width)
  
  p <- visNetwork(vis.nodes, vis.links, footer = cluster, ) %>%
    visIgraphLayout(layout = "layout_in_circle")
  
  return(p)
  
}

pathwayNetworkPlots <- function(ex.data, cls1, cls2=NULL, pathways, mgenes = NULL, geneType = "external_gene_name" ){
  library(manipulateWidget)
  library(htmltools)
  library(igraph)
  library(network)
  library(sna)
  library(ggraph)
  library(visNetwork)
  library(dplyr)
  library(qlcMatrix)
  library(data.table)
  library(stringr)
  library(readxl)
  library(KEGGREST)
  library(biomaRt)
  #ex.data <- dataFiltering(ex.data)
  if(!is.null(cls2)){
    expr.data <- ex.data[ex.data$phenotype %in% c(cls1, cls2), ]
  }else{
    expr.data <- ex.data
  }
  expr.data$phenotype <- ifelse(expr.data$phenotype == cls1, 0,1)
  cluster <- expr.data[expr.data$phenotype == 0,]
  other <- expr.data[expr.data$phenotype == 1,]
  
  print(c("Currently doing for pathway", pathways))
  genes = c()
  count = 0

  while(length(genes)==0 && count < 50){
    query <- tryCatch(keggGet(pathways), error = function(e){NULL})
    genes <- query[[1]]$GENE
    print(c(length(genes), count))
    count <- count+1
  }
  
  if(is.null(genes)){
    print(c("No genes identified for the pathway", pathways))
    next
  }
  print(c("KEGG Query done", pathways, geneType, organism = "mm"))
  genes <- genes[seq(2, length(genes), 2)]
  genes <- str_split(genes, ";", simplify = TRUE)[,1]
  genes <- toupper(genes)
  organism = "mm"
  if(organism == "mm"){
    ensmbl <- "mmusculus_gene_ensembl"
  }else if(organism == "hs"){
    ensmbl <- "hsapiens_gene_ensembl"
  }
  if(geneType == "ensembl_gene_id"){
    mart <- useEnsembl("ensembl", ensmbl)
    data_coords <-  getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filter = "external_gene_name",
                          values = genes, 
                          mart = mart)
    
    en_gen <- data_coords[data_coords$ensembl_gene_id%in% colnames(cluster),]$ensembl_gene_id
    ex_gen <- data_coords[data_coords$ensembl_gene_id%in% colnames(cluster),]$external_gene_name
  } else{
    
    en_gen <- genes[genes %in% colnames(ex.data)]
  }
   
  clsM <- as.matrix(subset(cluster, select = en_gen))
  otherM <- as.matrix(subset(other, select = en_gen))
    
  cls_res <- corSparse(clsM, clsM)
  oth_res <- corSparse(otherM, otherM)
  
  
  rownames(cls_res) <- colnames(clsM)
  colnames(cls_res) <- colnames(clsM)
  
  rownames(oth_res) <- colnames(otherM)
  colnames(oth_res) <- colnames(otherM)
  
  res_cls <- data.table("from" = rownames(cls_res), cls_res)
  res_cls <- melt(res_cls, id.vars = "from", value.name = "weights", variable.name = "to")
  
  res_cls[is.na(res_cls$weights)]$weights <- 0
  
  res_oth <- data.table("from" = rownames(oth_res), oth_res)
  res_oth <- melt(res_oth, id.vars = "from", value.name = "weights", variable.name = "to")
  
  res_oth[is.na(res_oth$weights)]$weights <- 0
  
  res_cls <- res_cls[res_cls$from != res_cls$to, ]
  res_oth <- res_oth[res_oth$from != res_oth$to, ]
  
  nodes <- data.frame("id" = en_gen)
  nodes$type <- ifelse(nodes$id %in% mgenes, "MARKER", "NORMAL")
  
  nodes_cls <- nodes
  nodes_oth <- nodes
  
  nodes_cls$cls_exp <- sapply(nodes$id, function(x) mean(clsM[,x]))
  nodes_oth$cls_exp <- sapply(nodes$id, function(x) mean(otherM[,x]))
  
  p1 <- plots(res_cls, nodes_cls,  cls1, "Cluster")
  p2 <- plots(res_oth, nodes_oth,  cls2, "Other")
  
  ff = combineWidgets(p1, p2, title = pathways)
  return(ff)
}


