library(Seurat)

# Creat Seurat Object
datList <- list()
folders <- c("WT1","WT2","KO1","KO2")
for(i in 1:length(folders))
  datList[[i]] <- Read10X(paste("/home/wll/projects/scTet20220707/02cellranger/",folders[i],"/outs/filtered_feature_bc_matrix/",sep=""))
sceList <- list()
for (i in 1:length(folders))
{
  tmp=CreateSeuratObject(counts = datList[[i]], project =folders[i],min.cells = 3, min.features = 200)
  tmp[["percent.mt"]]  <- PercentageFeatureSet(tmp, pattern = "^mt-")
  sceList[[i]] <- tmp
}
save(sceList,file="sceList.rda")