library(tidyverse)
library(patchwork)
library(Seurat)

#loading data
pbmc.data <- Read10X(data.dir = "C:/Users/tvb217/Documents/R/tmp/SCS/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#preprocessing
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
#normalizing data
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA<2500 &percent.mt < 5)
?subset

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


#identifying highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
plot2
#scaling data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#Linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#visualization in various modes
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = T)

#Determine dimensionality with null distribution
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

#alternative approach - Elbowplot
ElbowPlot(pbmc)

#clustering of cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc),5)

#non-linear dimensional reduction
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "C:/Users/tvb217/Documents/R/tmp/SCS/pbmc_tutorial.rds")

#finding markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers)
#finding markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0,3))
head(cluster5.markers,n = 5)

#finding markers for every cluster compared to all remaining cells, report only positives
pbmc.markers <- FindAllMarkers(pbmc, only.pos = T,
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
pbmc.markers %>% 
  group_by(cluster)%>%
  top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = T)
VlnPlot(pbmc, features = c("MS4A1", "CD79A", "NAMPT"))

#with row counts
VlnPlot(pbmc, features = c("NKG7", "PF4", "NAMPT"), slot = "counts", log = T)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

#Heatmap of expression for top20 markers
top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene)+ NoLegend()

#assigning cell type identity

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
?levels
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = T, pt.size = 0.5)+ NoLegend()

#saving the whole thing
saveRDS(pbmc, file = "C:/Users/tvb217/Documents/R/tmp/SCS/pbmc_final.rds")
