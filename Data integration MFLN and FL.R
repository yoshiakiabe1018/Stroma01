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
