#!/usr/bin/env Rscript
# Loading required packages
library(Seurat)
library(SeuratDisk)
library(anndata)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(dplyr)

# Set the variables
Input.dir = "Output/1_Label_Transfer_Annotation/"
Seurat.label = "predicted.celltype" #"Labelled_cluster"
Seurat.pattern = ".*Labelled_Clusters_SCT.1_UMAP_integrated_labelled.h5seurat"
seurat.list = list.files(Input.dir, pattern = Seurat.pattern, full.names = TRUE, recursive = TRUE)
print(seurat.list)

Predict.score.filtering <- function(seurat.x, score.cutoff = 0.5, annotation.type = "UMAP_integrated") {
  
  # Get the sample name
  seurat.sample = gsub(".*Output/.*//(.*?)/.*", "\\1", seurat.x)
  print(seurat.sample)
  
  # Set the output
  seurat.output = sub(".h5seurat", "_filtered.h5seurat", seurat.x)
  output.path = sub(".h5seurat", "", seurat.x)
  
  # Load the h5Seurat
  seurat.x = LoadH5Seurat(file = seurat.x)
  
  # Check how many bins passes the cutoff
  TF.table = table(seurat.x$predicted.celltype.score > score.cutoff)
  print(TF.table)
  
  # Generate a dataframe for downstream plotting
  predict.score.df = data.frame(seurat.x$predicted.celltype.score)
  predict.score.df$Sample = seurat.sample
  saveRDS(predict.score.df, file = paste0(output.path, "_predict_score.rds"))
  
  # Filter the seurat object
  seurat.x = subset(seurat.x, subset = predicted.celltype.score > score.cutoff)
  
  # Saved the filtered Seurat object
  SaveH5Seurat(object = seurat.x, filename = seurat.output, overwrite = TRUE)

  return(predict.score.df)
  
}

# Run the function
predict.score.filtering.out = lapply(seurat.list, Predict.score.filtering)
saveRDS(predict.score.filtering.out, file = "Output/1_Label_Transfer_Annotation/Predict_score_UMAP_Integrated_Major_Celltypes.rds")

