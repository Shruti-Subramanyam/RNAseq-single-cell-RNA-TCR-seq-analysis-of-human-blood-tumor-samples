#scRNA T cell w/doublet filtering 


MQcomb.markers %>%
  group_by(cluster) %>%library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)

MQpool1.data <- Read10X(data.dir = "E:/Bioinformatics/MQpool1_RS-03742653-20240117T054736Z-001/MQpool1_RS-03742653/count/sample_feature_bc_matrix")
MQpool2.data <- Read10X(data.dir = "E:/Bioinformatics/MQpool2_RS-03742653-20240117T054719Z-001/MQpool2_RS-03742653/count/sample_feature_bc_matrix")


MQpool1.obj <- CreateSeuratObject(counts = MQpool1.data$"Gene Expression",project = "MQpool1")
MQpool2.obj <- CreateSeuratObject(counts = MQpool2.data$"Gene Expression",project = "MQpool2")
#add.cell.ids specifies the new cell identity labels for the merged object
MQcomb <- merge(MQpool1.obj, y = MQpool2.obj, add.cell.ids = c("pool1","pool2"),project = "MQ_pool")

#combine different datatypes(layers are datatypes)
MQcomb <- JoinLayers(MQcomb)


head(colnames(MQcomb))
table(MQcomb$orig.ident)

#Calculating Percentage of Mitochondrial Genes
#PercentageFeatureSet is a Seurat function that calculates the percentage of mitochondrial (MT) genes for each cell
MQcomb[["percent.mt"]] <- PercentageFeatureSet(MQcomb, pattern = "^MT-")
#subset is used to filter cells based on specified criteria
MQcomb <- subset(MQcomb, nFeature_RNA >= 200 & nFeature_RNA <= 5000 & percent.mt < 15)

MQcomb <- NormalizeData(MQcomb)
#FindVariableFeatures identifies highly variable features (genes) in the dataset
#The function uses variance stabilizing transformation (VST) and selects the top 2000 variable features
MQcomb <- FindVariableFeatures(MQcomb, selection.method = "vst", nfeatures = 2000)#top 2000

all.genes <- rownames(MQcomb)
MQcomb <- ScaleData(MQcomb, features = all.genes) 
#RunPCA performs Principal Component Analysis on the dataset using the previously identified variable features
MQcomb <- RunPCA(MQcomb, features = VariableFeatures(object = MQcomb)) 

#RunUMAP computes the UMAP (dimensionality reduction) representation of the data using the first 20 principal components
MQcomb <- RunUMAP(MQcomb, dims = 1:20)
#FindNeighbors calculates a neighborhood graph based on the specified dimensions
MQcomb <- FindNeighbors(MQcomb, dims = 1:10)
#Find Clusters identifies cell clusters based on the neighborhood graph with a resolution of 0.5
MQcomb <- FindClusters(MQcomb, resolution = 0.5)
DimPlot(MQcomb, reduction = "umap", label = T)
#all the genes
all.genes <- rownames(MQcomb)
all.genes
#feature plot as per clusters
FeaturePlot(MQcomb, features = c("CD3E", "CD8A", "CD4", "FOXP3", "GATA3", "TBX21", "IL6ST",
                                 "TCF7", "SELL", "S1PR1", "CTLA4", "IL2RA", "IL7R", "CCR7",
                                 "GZMK", "CXCR3", "SLAMF6", "TNF", "IFNG", "CCL4", "CCL5",
                                 "PDCD1", "TIGIT", "HAVCR2", "GZMA", "GZMB", "PRF1", "CX3CR1"), pt.size = 0.5)
options(repr.plot.width = 8, repr.plot.height = 6)
FeaturePlot(MQcomb, features = c("CD3E", "CD8A", "CD4", "FOXP3", "GATA3", "TBX21", "IL6ST"), pt.size = 0.5)

#heatmap
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(MQcomb, features = top10$gene) + NoLegend()

#pK(proportion of doublets) Identification
#identify and remove doublets from single-cell RNA sequencing
#paramSweep is a function from doubletFinder that performs a parameter sweep across different values of the pK parameter, which represents the proportion of doublets in the dataset
sweep.res.list_MQ <- paramSweep(MQcomb, PCs = 1:20, sct = FALSE)
#GT = FALSE indicates that the ground truth for doublets is not provided
sweep.state_MQ <- summarizeSweep(sweep.res.list_MQ, GT = FALSE)

#find.pK is used to identify the optimal pK value based on the summarized sweep results
bcmvn_MQ.combined <- find.pK(sweep.state_MQ)
ggplot(bcmvn_MQ.combined, aes(pK, BCmetric, group = 1)) +geom_point() + geom_line() + RotatedAxis()

#pk.set <- unique(sweep.state_MQ$pK)[2]


nExp_poi <- round(0.08*nrow(MQcomb@meta.data))
#MQcomb <- doubletFinder(MQcomb, PCs = 1:20, pN = 0.25, pK = as.numeric(as.character(pk.set)),nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#doubletFinder is used to detect doublets in the scRNA-seq data
#The function takes parameters such as the top PCs, pN (expected proportion of doublets), pK (proportion of doublets in the negative control), nExp (number of expected doublets), and other options
MQcomb <- doubletFinder(MQcomb, PCs = 1:20, pN = 0.25, pK = 0.28,nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(MQcomb)

#MQcomb <- doubletFinder(MQcomb, PCs = 1:20, pN = 0.25, pK = 0.01,nExp = nExp_poi, reuse.pANN = "paNN_0.25_0.01_272", sct = FALSE)
#Plot for Doublets??DF.classifications_0.25_0.28_626??# set it as the last column of the MQcomb table
a <- DimPlot(MQcomb,pt.size = 1,label=TRUE, label.size = 5,reduction = "umap",group.by = "DF.classifications_0.25_0.28_626" )+theme(aspect.ratio = 1)
a
#Remaining single cells
MQcomb_doub_rem <- subset(MQcomb, subset = DF.classifications_0.25_0.28_626 == "Singlet")
head(MQcomb_doub_rem)

#number of singlets and doublets
table(MQcomb@meta.data$DF.classifications_0.25_0.28_626)

write_csv(MQcomb@meta.data,"E:/Bioinformatics/MQcomb_data.csv")

#analysis markers

#Normalization and Variable Feature Selection
MQcomb_doub_rem <- NormalizeData(MQcomb_doub_rem)
MQcomb_doub_fea <- FindVariableFeatures(MQcomb_doub_rem, selection.method = "vst", nfeatures = 2000)#variance stabilizing transformation
#Variable Feature Plot with Top Features
top10 <- head(VariableFeatures(MQcomb_doub_fea), 10)
plot1 <- VariableFeaturePlot(MQcomb_doub_fea)
#Error??
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

#Scaling Data and Running PCA
all.genes.ana <- rownames(MQcomb_doub_fea)
MQcomb_doub_fea <- ScaleData(MQcomb_doub_fea, features = all.genes.ana)

MQcomb_doub_fea <- RunPCA(MQcomb_doub_fea, features = VariableFeatures(object = MQcomb_doub_fea))
print(MQcomb_doub_fea[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(MQcomb_doub_fea, dims =4, reduction = "pca")
DimPlot(MQcomb_doub_fea, reduction = "pca") + NoLegend()
DimHeatmap(MQcomb_doub_fea, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(MQcomb_doub_fea, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(MQcomb_doub_fea)

MQcomb_doub_fea <- FindNeighbors(MQcomb_doub_fea, dims = 1:10)
MQcomb_doub_fea <-FindClusters(MQcomb_doub_fea, resolution = 0.5)

MQcomb_doub_fea <-RunUMAP(MQcomb_doub_fea, dims = 1:10)
DimPlot(MQcomb_doub_fea, reduction = "umap", label = TRUE)

cluster2.markers <- FindMarkers(MQcomb_doub_fea, ident.1 = 2)
head(cluster2.markers, n =5)

#CDI sample gene - as there is no variation, so no graph was displayed.
VlnPlot(MQcomb_doub_fea, features = c("S100A11","IL32"))
FeaturePlot(MQcomb_doub_fea, features = c('CDY1'))
