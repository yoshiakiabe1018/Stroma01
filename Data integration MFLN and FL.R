library(dplyr)
library(Seurat)
library(patchwork)


=== DATA INTEGRATION ===

anchors_MFLN9_20 <- FindIntegrationAnchors(object.list = list(MFLN01,
                                                              MFLN02,
                                                              MFLN03,
                                                              MFLN04,
                                                              MFLN05,
                                                              MFLN06,
                                                              MFLN07,
                                                              MFLN08,
                                                              MFLN09),
                                             dims = 1:20)
Combined_MFLN9_20 <- IntegrateData(anchorset = anchors_MFLN9_20, dims = 1:20)

anchors_FL10_20 <- FindIntegrationAnchors(object.list = list(FL01,
                                                             FL02,
                                                             FL03,
                                                             FL04,
                                                             FL05,
                                                             FL06,
                                                             FL07,
                                                             FL08,
                                                             FL09,
                                                             FL10),
                                          dims = 1:20)
Combined_FL10_20 <- IntegrateData(anchorset = anchors_FL10_20, dims = 1:20)

anchors_9vs10_20 <- FindIntegrationAnchors(object.list = list(Combined_MFLN9_20,
                                                              Combined_FL10_20),
                                                           dims = 1:20)
Combined_9vs10_20 <- IntegrateData(anchorset = anchors_9vs10_20, dims = 1:20)


=== SCALING AND PCA ===

Combined_9vs10_20_scaled <- ScaleData(Combined_9vs10_20, verbose = FALSE)
Combined_9vs10_20_PCA <- RunPCA(Combined_9vs10_20_scaled, npcs = 200, verbose = FALSE)


=== CLUSTERING ===

Combined_9vs10_20_PCA_180 <- RunUMAP(Combined_9vs10_20_PCA, reduction = "pca", dims = 1:180)
Combined_9vs10_20_PCA_180 <- FindNeighbors(Combined_9vs10_20_PCA_180, reduction = "pca", dims = 1:180)
Combined_9vs10_20_PCA_180_1.8 <- FindClusters(Combined_9vs10_20_PCA_180, resolution = 1.8)


#After supervised annotation
MAJOR3 <- subset(MAJOR6,
                 idents = c("BEC","LEC","NESC"),
                 invert = F)  

MAJOR3_MFLN <- subset(x = MAJOR3, subset = FLvsMFLN == "MFLN")

MAJOR3_MFLN_averages <- AverageExpression(MAJOR3_MFLN)
orig_levels_MAJOR3_MFLN_averages <- levels(MAJOR3_MFLN)
Idents(MAJOR3_MFLN) <- gsub(pattern = " ", replacement = "_", x = Idents(MAJOR3_MFLN))
orig_levels_MAJOR3_MFLN_averages <- gsub(pattern = " ", replacement = "_", x = orig_levels_MAJOR3_MFLN_averages)
levels(MAJOR3_MFLN) <- orig_levels_MAJOR3_MFLN_averages
Cluster_averages_MAJOR3_MFLN <- AverageExpression(MAJOR3_MFLN, return.seurat = TRUE)

DoHeatmap(Cluster_averages_MAJOR3_MFLN,
          features = Top30DEGs_MAJOR3_MFLN,
          size = 4.5, 
          disp.min = -2.5,
          disp.max = 3,
          slot = "scale.data",
          angle = 0,
          draw.lines = FALSE,
          group.bar.height = 0) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))

  
=== SUBCLUSTERING ANALYSIS ===

BEC <- subset(MAJOR3,idents = c("BEC"),invert = F)
LEC <- subset(MAJOR3,idents = c("LEC"),invert = F)
NESC <- subset(MAJOR3,idents = c("NESC"),invert = F)
DefaultAssay(BEC) <- "integrated"
DefaultAssay(LEC) <- "integrated"
DefaultAssay(NESC) <- "integrated"
BEC_scaled <- ScaleData(BEC, verbose = FALSE)
LEC_scaled <- ScaleData(LEC, verbose = FALSE)
NESC_scaled <- ScaleData(NESC, verbose = FALSE)
BEC_PCA <- RunPCA(BEC_scaled, npcs = 100, verbose = FALSE)
LEC_PCA <- RunPCA(LEC_scaled, npcs = 100, verbose = FALSE)
NESC_PCA <- RunPCA(NESC_scaled, npcs = 100, verbose = FALSE)

BEC_PCA_60 <- RunUMAP(BEC_PCA, reduction = "pca", dims = 1:60)
BEC_PCA_60 <- FindNeighbors(BEC_PCA_60, reduction = "pca", dims = 1:60)
BEC_PCA_60_0.5 <- FindClusters(BEC_PCA_60, resolution = 0.5)

LEC_PCA_80 <- RunUMAP(LEC_PCA, reduction = "pca", dims = 1:80)
LEC_PCA_80 <- FindNeighbors(LEC_PCA_80, reduction = "pca", dims = 1:80)
LEC_PCA_80_0.5 <- FindClusters(LEC_PCA_80, resolution = 0.5)

NESC_PCA_60 <- RunUMAP(NESC_PCA, reduction = "pca", dims = 1:60)
NESC_PCA_60 <- FindNeighbors(NESC_PCA_60, reduction = "pca", dims = 1:60)
NESC_PCA_60_1.0 <- FindClusters(NESC_PCA_60, resolution = 1.0)

#QC
