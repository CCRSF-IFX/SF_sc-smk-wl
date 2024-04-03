library(Seurat)
library(Matrix)
library(MASS)
library(dplyr)
library(Biostrings)
library(reshape2)

seurat_obj <- CreateSeuratObject(counts = Read10X(data.dir = path2cnt_gene_BR1_1), project = "MAS_BR1_1")

