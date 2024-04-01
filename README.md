# scSCOPE
It is challenging to select marker genes robustly to discriminate cell types with meaningful cellular functions. Here, we propose an extension of our bulk RNA-seq tool, Stabilized COre gene and Pathway Election, or SCOPE, to identify marker genes and pathways in single-cell RNA-seq data (scRNA-seq). The new tool, single-cell-SCOPE, or scSCOPE, integrates LASSO, Co-expression analyses, Pathway Enrichment and Differential Expression to identify marker genes and pathways in single-cell data. 

![Fig1_new](https://github.com/QingrunZhangLab/scSCOPE/assets/66895308/ea6e3df0-abf1-45c4-b7e4-c9df31468e62)

The input for scSCOPE is a gene expression matrix with cells in rows, genes in columns and phenotype/cluster information for each cell in a column titled "phenotype". Once you have saved the file in csv format, you can follow scSCOPE_WorkFlow to identify marker genes and pathways in your data. 
