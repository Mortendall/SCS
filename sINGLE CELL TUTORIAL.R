library(tidyverse)
library(patchwork)
library(Seurat)

pbmc.data <- Read10X(data.dir = "C:/Users/tvb217/Documents/R/SCS/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
