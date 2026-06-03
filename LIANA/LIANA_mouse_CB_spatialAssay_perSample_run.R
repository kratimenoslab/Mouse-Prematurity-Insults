library(spatstat.explore)
library(Seurat)
library(liana)
library(dplyr)
library(stringr)
library(ggplot2)

## LOAD DATA
load('mouse_msobj_transformed_analysis.rda')
outputDir <- 'liana/'
prefix <- 'mouseSpatial_PCGC_perSample_'


#####################################################
################# PREPPING DATA #####################
#####################################################

subsur <- msobj 
DefaultAssay(subsur) <- 'Spatial' # set spatial assay, [1] 31053  8961
# checked counts and data are the same

# Subset spots to PC GC only
subsur <- subsur[, subsur$annotations %in% c("PC", "GC")] # subset to PC GC only, [1] 31053  2430
print(dim(subsur)) # [1] 31053  2430

# Remove genes with no signal
counts_mat <- GetAssayData(subsur, assay = "Spatial", layer = "counts")
keep_genes <- apply(counts_mat, 1, max) > 0
subsur <- subsur[keep_genes, ] 
print(dim(subsur)) # [1] 19792  2430

### ortho vignette https://github.com/saezlab/liana/blob/6cab46c54234f861ea176c3de77c4b8aa45ecb3d/vignettes/liana_ortho.Rmd#L33
# Here, we will convert LIANA's Consensus resource to murine symbols
op_resource <- select_resource("Consensus")[[1]]

# Generate orthologous resource
ortholog_resource <- generate_homologs(op_resource = op_resource,
                                       target_organism = 10090) # mouse


#####################################################
################# RUNNING LIANA #####################
#####################################################
alllianas_agg <- data.frame()

for (pan in unique(subsur$uniquepanels)){
  
  print(pan)
  testdata <- subsur[,subsur$uniquepanels == pan] # subset to the panel

  # Check both cell types are present with enough cells
  cell_counts <- table(testdata$annotations)
  if (!all(c("PC", "GC") %in% names(cell_counts)) ||
      any(cell_counts[c("PC", "GC")] < 5)) {
    message("Skipping ", pan, " — insufficient cells")
    next
  }

  testdata <- NormalizeData(testdata) # this normalizes the 'data' slot, not the 'counts' slot; log1p(counts / total_counts_per_cell * scale.factor)
  Idents(testdata) <- "annotations"

  liana_test <- liana_wrap(testdata, assay="Spatial", min_cells = 5, resource = "custom", idents_col='annotations',
                           external_resource = ortholog_resource)

  liana_test <- liana_test %>% liana_aggregate() #aggregates results for each slide
  liana_test$uniquepanel <- pan
  alllianas_agg <- rbind(alllianas_agg, as.data.frame(liana_test))

}

saveRDS(alllianas_agg, paste0(outputDir, prefix, 'alllianasagg.RDS')) #aggregated per slide
write.csv(alllianas_agg, paste0(outputDir, prefix, 'alllianasagg.csv'))


