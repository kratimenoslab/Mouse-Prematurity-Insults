library(CoGAPS)
library(Matrix)
library(gdata)
library(Seurat)
library(ggplot2)

## Run for all spots
load("/mnt/morbo/Data/Users/kwoyshner/cerebellum/data/mouse_msobj_transformed_analysis.rda")
spatial_matrix <- as.matrix(msobj@assays[["Spatial"]]@counts) # genes by spots
spatial_assay <- msobj@assays[['Spatial']] # genes x cells

# Remove spots with no signal
spatial_matrix <- spatial_matrix[,apply(spatial_matrix,2,max)>0]
# Remove genes with no signal
spatial_matrix <- spatial_matrix[apply(spatial_matrix,1,max)>0,]

log_spatial_matrix <- log1p(spatial_matrix)

## Set params
nPatterns = 15
nIterations = 20000

params <- new("CogapsParams")
geneNames <- rownames(log_spatial_matrix)
spotNames <- colnames(log_spatial_matrix)
params <- CogapsParams(
  sparseOptimization=FALSE,
  nPatterns=nPatterns,
  seed=123,
  geneNames=geneNames,
  sampleNames=spotNames,
  nIterations=20000,
  distributed='genome-wide'
)

cogaps.exprs<-log_spatial_matrix
params <- setDistributedParams(params, nSets=7)

outputDir <- '/mnt/morbo/Data/Users/kwoyshner/cerebellum/results/mouse_allGenes_n15/'
savename <- paste0('mouse_CB_cogaps_n',nPatterns,'_nIterations',nIterations/1000,'k')


NMF<-CoGAPS(as.matrix(cogaps.exprs),params = params, distributed="genome-wide", outputFrequency=500)

saveRDS(NMF,file=paste0(outputDir,savename,".RDS"))

sampleWeights<-t(NMF@sampleFactors)
colnames(sampleWeights)<-colnames(geneExp)
write.csv(sampleWeights, paste0(outputDir, savename, "_patterns.csv"))

geneWeights<-NMF@featureLoadings
rownames(geneWeights)<-rownames(geneExp)
write.csv(geneWeights, paste0(outputDir, savename, "_geneWeights.csv"))

