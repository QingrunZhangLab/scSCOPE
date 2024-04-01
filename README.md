# scSCOPE
It is challenging to select marker genes robustly to discriminate cell types with meaningful cellular functions. Here, we propose an extension of our bulk RNA-seq tool, Stabilized COre gene and Pathway Election, or SCOPE, to identify marker genes and pathways in single-cell RNA-seq data (scRNA-seq). The new tool, single-cell-SCOPE, or scSCOPE, integrates LASSO, Co-expression analyses, Pathway Enrichment and Differential Expression to identify marker genes and pathways in single-cell data. 

![Fig1_new-min](https://github.com/QingrunZhangLab/scSCOPE/assets/66895308/0503b5be-013a-4e85-8f9d-f06e48956b0f)

The input for scSCOPE is a gene expression matrix with cells in rows, genes in columns and phenotype/cluster information for each cell in a column titled "phenotype". Once you have saved the file in csv format, you can follow scSCOPE_WorkFlow to identify marker genes and pathways in your data.

For example: The code below can find marker genes between cluster 1 and cluster 5. We recommend to set the "iter" and "iter_lasso" to default value. 
Detailed explanation of each functions and parameters are given below. 
```
findMarkerGenes(ex.data = ex.data, cluster2 = 5, 
                cluster1 = 1, iters = 10, iters_lasso = 10, org = "mmusculus", pathDB = "KEGG", 
                iter_correlation = 10, core_cutoff = 8, FC_threshold = 0.5, corr_iter_cut = 8, 
                geneId = "ensembl_gene_id")
```

