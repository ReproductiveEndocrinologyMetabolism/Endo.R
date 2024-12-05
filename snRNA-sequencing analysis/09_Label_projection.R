#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(cowplot)

# Generate or set the output directory
if (dir.exists(path = "Output/9_Main_relabelled") == FALSE) {
  print("Output/9_Main_relabelled")
  dir.create(path = "Output/9_Main_relabelled", recursive = TRUE)
  Output.dir = "Output/9_Main_relabelled/"
} else if (dir.exists(path = "Output/9_Main_relabelled") == TRUE) {
  print("Directory exists")
  Output.dir = "Output/9_Main_relabelled/"
} else {
  print("Error with output directory")
}

Project_name = "Endo_All"
celltype.pattern = ".*Endo_All_(.*?)_reclustered_labelled\\.h5seurat"
seurat.main.path = "Output/3_Labeling/Endo_All_SCT_celltypes_integrated.h5seurat"

# List the labelled seurat objects to be projected on the main
seurat.dirs = list.dirs(path = "Output")
seurat.dirs = seurat.dirs[grep("Output/5_.*_Clustering", seurat.dirs)]
print(seurat.dirs)

# For-loop to extract labels from each labelled cell type
seurat.files = list.files(seurat.dirs, pattern = "Endo_All_.*_reclustered_labelled.h5seurat", full.names = TRUE)
print(seurat.files)

Extract.seurat.label <- function(seurat.x) {
  
  # Extract which celltype that is processed
  celltype.x = gsub(celltype.pattern, "\\1", seurat.x)
  print(paste("Extraxting label from", celltype.x))
  
  # Load the seurat object
  seurat.x = LoadH5Seurat(file = seurat.x)
  
  # Extract the label
  seurat.label = seurat.x[[paste0(celltype.x, "_labelled")]]
  colnames(seurat.label) = "Label"
  
  # Save the celltype labels
  saveRDS(seurat.label, file = paste0(Output.dir, Project_name, "_", celltype.x, "_Labels.rds"))
  
  # Return the label
  return(seurat.label)
  
}

# Run the function
labels.list = lapply(seurat.files, Extract.seurat.label)

# Combine the labels to one combined dataframe and save it
labels.combined = do.call(rbind, labels.list)
saveRDS(labels.combined, file = paste0(Output.dir, Project_name, "_Labels_combined.rds"))

# Load the original main seurat object and label it.
seurat.main = LoadH5Seurat(file = seurat.main.path)
seurat.main$Combined_labels = labels.combined
SaveH5Seurat(seurat.main, filename = paste0(Output.dir, Project_name, "_Combined_labels.h5seurat"), overwrite = TRUE)

# Visualise the results by UMAP
UMAP.x = DimPlot(seurat.main, reduction = "umap", group.by = "Combined_labels",
                 label = TRUE, repel = TRUE)
ggsave2(file = paste0(Output.dir, Project_name, "_UMAP_LABELLED_Combined_labels.pdf"), plot = UMAP.x, dpi = 700)

UMAP.x = DimPlot(seurat.main, reduction = "umap", group.by = "Combined_labels")
ggsave2(file = paste0(Output.dir, Project_name, "_UMAP_UNLABELLED_Combined_labels.pdf"), plot = UMAP.x, dpi = 700)

print("DONE")
