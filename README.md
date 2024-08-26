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


Gene Network Plots for a marker gene of interest can be created by using geneNetworkPlots() function.
```
all_p <- geneNetworkPlots(expr.data, marker_gene, cluster1, cluster2= NULL, pathways_of_interest)![Uploading Screenshot 2024-08-20 at 7.48.26 PM.png…]()

#For example: 
all_p <- geneNetworkPlots(expr.data, "Il7r", cluster1 = 4, cluster2 = 10, pathways_of_interest = c("mmu04660", "mmu04659"))
print(all_p$Il7r)
```
<geneNetwork src="https://github.com/user-attachments/assets/80c9fc4c-4f1c-4882-8af7-b2698c2979b5" width="100" height="100")>


Pathway Network Plot for a pathway of interest in two clusters can be created by using pathwayNetworkPlots() function. 

```
pathwayNetworkPlot(ex.data, cls1, cls2, pathways, geneType = "ensembl_gene_id")

#For example:
p <- pathwayNetworkPlots(ex.data, cls1 = 5, cls2 = 4, pathways = "mmu04659", geneType = "ensembl_gene_id")
print(p)

```

<img width="201" alt="Screenshot 2024-08-20 at 7 48 26 PM" src="https://github.com/user-attachments/assets/0dc4cc66-59eb-48a6-9d75-79206bd61d85">

