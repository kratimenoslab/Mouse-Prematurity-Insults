library(CoGAPS)
library(Matrix)
library(gdata)
library(Seurat)
library(ggplot2)

## Run for all spots subset to folia
load("/mnt/morbo/Data/Users/kwoyshner/cerebellum/data/mouse_msobj_transformed_analysis.rda")
msobj_folia <- msobj[,msobj$folia == "folia"]
print(dim(msobj_folia)) # [1] 19372  5830;

spatial_matrix <- as.matrix(msobj_folia@assays[["Spatial"]]@counts) # genes by spots
print(dim(spatial_matrix)) # [1] 31053  8961; genes x spots

# Remove spots with no signal
spatial_matrix <- spatial_matrix[,apply(spatial_matrix,2,max)>0]
# Remove genes with no signal
spatial_matrix <- spatial_matrix[apply(spatial_matrix,1,max)>0,]
print(dim(spatial_matrix))

log_spatial_matrix <- log1p(spatial_matrix)

## Set params
nPatterns = 15
nIterations = 20000
outputDir <- 'mouse_cogaps_updated/'
savename <- paste0('mouse_CB_folia_cogaps_n',nPatterns,'_nIterations',nIterations/1000,'k')

print(paste0('Saving output to ', outputDir, savename,".RDS")) # This is the file we're saving to
cat("Test",file=paste0(outputDir,savename,"_test.txt"),append=TRUE) # Save a test file here to make sure we can write to this directory

params <- new("CogapsParams")
geneNames <- rownames(log_spatial_matrix)
spotNames <- colnames(log_spatial_matrix)
params <- CogapsParams(
  sparseOptimization=FALSE,
  nPatterns=nPatterns,
  seed=123,
  geneNames=geneNames,
  sampleNames=spotNames,
  nIterations=nIterations,
  distributed='genome-wide'
)

cogaps.exprs<-log_spatial_matrix
params <- setDistributedParams(params, nSets=5)

NMF<-CoGAPS(as.matrix(cogaps.exprs),params = params, distributed="genome-wide", outputFrequency=1000)
saveRDS(NMF,file=paste0(outputDir,savename,".RDS"))

sampleWeights<-NMF@sampleFactors
#colnames(sampleWeights)<-colnames(geneExp)
write.csv(sampleWeights, paste0(outputDir, savename, "_patterns.csv"))

geneWeights<-NMF@featureLoadings
#rownames(geneWeights)<-rownames(geneExp)
write.csv(geneWeights, paste0(outputDir, savename, "_geneWeights.csv"))

