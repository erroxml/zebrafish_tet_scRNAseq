library(dplyr)
library(Seurat)

load("scTet_integrate.rda")
load("sceList.rda")
folders <- c("WT1","WT2","KO1","KO2")
#select NC cells and reanalysis
ident_df <- data.frame(cell=names(Idents(scTet_integrate)), cluster=Idents(scTet_integrate),sid=scTet_integrate@meta.data$orig.ident)

sce.norm.sub <- list()
for(i in 1:length(folders))
{
  sce.norm.sub[[i]] <- RenameCells(sceList[[i]], add.cell.id = folders[i])
  sce.norm.sub[[i]] <- subset(sce.norm.sub[[i]], subset = nFeature_RNA>1200 & nFeature_RNA<5000 & percent.mt<5  & nCount_RNA<20000)
  sce.norm.sub[[i]] <- subset(sce.norm.sub[[i]], cells=as.vector(ident_df[which(ident_df$cluster%in% c(4,7,14)),1]))
  sce.norm.sub[[i]] <- NormalizeData(sce.norm.sub[[i]], verbose = FALSE)
  sce.norm.sub[[i]] <- FindVariableFeatures(sce.norm.sub[[i]], selection.method = "vst", nfeatures = 2000)
}
anchors.sub <- FindIntegrationAnchors(object.list = sce.norm.sub)
scTet_integrate_filter <- IntegrateData(anchorset = anchors.sub)
scTet_integrate_filter <- ScaleData(scTet_integrate_filter, verbose = FALSE)
scTet_integrate_filter <- RunPCA(scTet_integrate_filter, npcs = 50, verbose = FALSE)
scTet_integrate_filter <- FindNeighbors(scTet_integrate_filter, reduction = "pca", dims = 1:29)
scTet_integrate_filter <- FindClusters(scTet_integrate_filter, resolution = .1)
scTet_integrate_filter <- RunUMAP(scTet_integrate_filter, reduction = "pca", dims = 1:29)
scTet_integrate_filter <- RunTSNE(scTet_integrate_filter, reduction = "pca", dims = 1:29)
scTet_integrate_filter <- FindVariableFeatures(scTet_integrate_filter, selection.method = "vst", nfeatures = 2000)

save(scTet_integrate_filter,file = "scTet_integrate_filter.rda")
sub_markers <- FindAllMarkers(object = scTet_integrate_filter, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.2)
write.csv(sub_markers,file = "sub_marker.csv")

#unselect cluster4(neurons) and reanalysis
ident_df2 <- data.frame(cell=names(Idents(scTet_integrate_filter)), cluster=Idents(scTet_integrate_filter),sid=scTet_integrate_filter@meta.data$orig.ident)

sce.norm.sub2 <- list()
for(i in 1:length(folders))
{
  sce.norm.sub2[[i]] <- RenameCells(sceList[[i]], add.cell.id = folders[i])
  sce.norm.sub2[[i]] <- subset(sce.norm.sub2[[i]], subset = nFeature_RNA>1200 & nFeature_RNA<5000 & percent.mt<5  & nCount_RNA<20000)
  sce.norm.sub2[[i]] <- subset(sce.norm.sub2[[i]], cells=as.vector(ident_df2[which(ident_df2$cluster%in% c(0,1,2,3,5,6,7,8)),1]))
  sce.norm.sub2[[i]] <- NormalizeData(sce.norm.sub2[[i]], verbose = FALSE)
  sce.norm.sub2[[i]] <- FindVariableFeatures(sce.norm.sub2[[i]], selection.method = "vst", nfeatures = 2000)
}
anchors.sub2 <- FindIntegrationAnchors(object.list = sce.norm.sub2)
scTet_integrate_filter2 <- IntegrateData(anchorset = anchors.sub2)
scTet_integrate_filter2 <- ScaleData(scTet_integrate_filter2, verbose = FALSE)
scTet_integrate_filter2 <- RunPCA(scTet_integrate_filter2, npcs = 50, verbose = FALSE)
scTet_integrate_filter2 <- FindNeighbors(scTet_integrate_filter2, reduction = "pca", dims = 1:22)
scTet_integrate_filter2 <- FindClusters(scTet_integrate_filter2, resolution = 1)
scTet_integrate_filter2 <- RunUMAP(scTet_integrate_filter2, reduction = "pca", dims = 1:22)
scTet_integrate_filter2 <- RunTSNE(scTet_integrate_filter2, reduction = "pca", dims = 1:22)
scTet_integrate_filter2 <- FindVariableFeatures(scTet_integrate_filter2, selection.method = "vst", nfeatures = 2000)

save(scTet_integrate_filter2,file = "scTet_intergate_filter2.rda")
sub_markers2 <- FindAllMarkers(object = scTet_integrate_filter2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(sub_markers2,file = "sub_marker2.csv")

#unselect cluster6,7,8,9(epilsdermis or others)
scTet_integrate_filter.onlync <- subset(scTet_integrate_filter2,subset=seurat_clusters %in% c(0,1,2,3,4,5))
scTet_integrate_filter.onlync.markers <- sub_markers2[sub_markers2$cluster %in% c(0,1,2,3,4,5),]
save(scTet_integrate_filter.onlync, file = "scTet_integrate_filter.onlync.rda")
write.csv(scTet_integrate_filter.onlync.markers, file = "scTet_integrate_filter.onlync.markers.csv")