library(DoubletFinder)
library(Seurat)
library(dplyr)
load("scTet_integrate.rda")

#Find the optimal pK value
sweep.res.list <- DoubletFinder::paramSweep_v3(seu = scTet_integrate, PCs = 1:23, sct = TRUE)
sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
pk_scTet <- DoubletFinder::find.pK(sweep.stats = sweep.stats)
pk <- pk_scTet$pK[which.max(pk_scTet$BCmetric)] %>% as.character() %>% as.numeric()
## Optimize the expected number of doublets
DoubletRate = 0.06 #~8000 cells
homotypic.prop <- DoubletFinder::modelHomotypic(annotations = scTet_integrate@meta.data$annotation)          
nExp_poi <- round(DoubletRate*nrow(scTet_integrate@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

scTet_integrate_doublet_output <- DoubletFinder::doubletFinder_v3(seu = scTet_integrate, 
                                                        PCs = 1:23, 
                                                        pK = pk, 
                                                        nExp = nExp_poi.adj, 
                                                        reuse.pANN = FALSE,
                                                        sct = TRUE)
save(scTet_integrate_doublet_output,file = "scTet_integrate_doublet_output.rda")

DimPlot(scTet_integrate_doublet_output, group.by = "DF.classifications_0.25_0.09_1234", reduction = "umap") + 
  DimPlot(scTet_integrate, group.by = "annotation", reduction = "umap")

FeaturePlot(scTet_integrate_doublet_output, "pANN_0.25_0.09_1234", reduction = "umap")+ 
  DimPlot(scTet_integrate, group.by = "annotation", reduction = "umap")

VlnPlot(scTet_integrate_doublet_output, features = "pANN_0.25_0.09_1234", group.by = "annotation")