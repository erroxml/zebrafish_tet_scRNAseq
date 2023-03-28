library(dplyr)
library(Seurat)
library(ggplot2)

load("sceList.rda")
# filter and normalize
folders <- c("WT1","WT2","KO1","KO2")
sce.norm <- list()
for(i in 1:length(folders))
{
  sce.norm[[i]] <- RenameCells(sceList[[i]], add.cell.id = folders[i])
  sce.norm[[i]] <- subset(sce.norm[[i]], subset = nFeature_RNA >1200 & nFeature_RNA<5000 & percent.mt<5  & nCount_RNA<20000)
  sce.norm[[i]] <- NormalizeData(sce.norm[[i]], verbose = FALSE)
  sce.norm[[i]] <- FindVariableFeatures(sce.norm[[i]], selection.method = "vst", nfeatures = 2000)
}

#integrate
anchors <- FindIntegrationAnchors(object.list = sce.norm)
scTet_integrate <- IntegrateData(anchorset = anchors)

#DefaultAssay(scTet_integrate) <- "integrated"

# Run the standard workflow for visualization and clustering
scTet_integrate <- ScaleData(scTet_integrate, verbose = FALSE)
scTet_integrate <- RunPCA(scTet_integrate, npcs = 50, verbose = FALSE)

#ElbowPlot(scTet_integrate,ndims = 50)
#t-SNE and Clustering
scTet_integrate <- FindNeighbors(scTet_integrate, reduction = "pca", dims = 1:23)
scTet_integrate <- FindClusters(scTet_integrate, resolution = .2)
scTet_integrate <- RunUMAP(scTet_integrate, reduction = "pca", dims = 1:23)
scTet_integrate <- RunTSNE(scTet_integrate, reduction = "pca", dims = 1:23)
scTet_integrate <- FindVariableFeatures(scTet_integrate, selection.method = "vst", nfeatures = 2000)

###find marker genes for each cluster
markers <- FindAllMarkers(object = scTet_integrate, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(markers,file = "marker.csv")

#annotation
metadata <- scTet_integrate@meta.data
metadata$type <- substr(metadata$orig.ident,start = 1L,stop = 2L)
metadata$annotation <- "Undefined"
metadata$annotation[which(metadata$seurat_clusters %in% c(0,3))] <- "Nerual Progenitor cell"
#metadata$annotation[which(metadata$seurat_clusters %in% c(3))] <- "proliferative region"
metadata$annotation[which(metadata$seurat_clusters %in% c(1,2,6,9))] <- "Differentiating neuron"
metadata$annotation[which(metadata$seurat_clusters %in% c(5))] <- "Eye Photorecepter cell"
metadata$annotation[which(metadata$seurat_clusters %in% c(4,7,14))] <- "Ectomesenchymal"
metadata$annotation[which(metadata$seurat_clusters %in% c(8,17))] <- "Epidermal"
metadata$annotation[which(metadata$seurat_clusters %in% c(10))] <- "Erythrocyte"
metadata$annotation[which(metadata$seurat_clusters %in% c(11))] <- "Muscle"
metadata$annotation[which(metadata$seurat_clusters %in% c(12))] <- "Xanthophore"
metadata$annotation[which(metadata$seurat_clusters %in% c(13))] <- "Melanocyte"
metadata$annotation[which(metadata$seurat_clusters %in% c(15))] <- "Hemangioblast"
metadata$annotation[which(metadata$seurat_clusters %in% c(18))] <- "Iridophore"
metadata$annotation[which(metadata$seurat_clusters %in% c(16,19))] <- "Phagocytes"
scTet_integrate@meta.data <- metadata
save(scTet_integrate,file = "scTet_integrate.rda")

annotation_markers <- data.frame()
for(i in 1:length(unique(metadata$annotation))){
  tmp_markers <- FindMarkers(object = scTet_integrate, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2,group.by="annotation",ident.1=unique(metadata$annotation)[i])
  tmp_markers$cluster <- unique(metadata$annotation)[i]
  annotation_markers <- rbind(annotation_markers,tmp_markers)
}
annotation_markers$gene <- rownames(annotation_markers)
write.csv(annotation_markers,file = "annotation_markers.csv")
top20 <- annotation_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(scTet_integrate, features = top20$gene) + NoLegend()
