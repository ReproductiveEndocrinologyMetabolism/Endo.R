#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(stringr)
library(openxlsx)
library(pheatmap)

#celltype = "Stromal"
celltype = "Immune" # Epithelium, Immune, Endothelial
Project_name = paste0("Endo_All_", celltype)
#Project_name = "Endo_All_Stromal_uSMC"
Input.dir = paste0("Output/5_", celltype, "_Clustering/")
#Input.dir = paste0("Output/5_Stromal_uSMC_Clustering/")
Norm.assay = "RNA"
CtrlvsPCOS.flag = TRUE

if (dir.exists(path = paste0("Output/8_DEG_GO_plotting_", celltype)) == FALSE) {
  print(paste0("Generating output directory Output/8_DEG_GO_plotting_", celltype))
  dir.create(path = paste0("Output/8_DEG_GO_plotting_", celltype), recursive = TRUE)
  Output.dir = paste0("Output/8_DEG_GO_plotting_", celltype, "/")
} else if (dir.exists(path = paste0("Output/8_DEG_GO_plotting_", celltype)) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/8_DEG_GO_plotting_", celltype, "/")
} else {
  print("Error with output directory")
}

# Loading reclustered Seurat
print("Seurat object loading")
x.seurat = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_reclustered_labelled.h5seurat"))

# If CtrlvsPCOS.flag is true, only baseline will be subset for further scaling and plotting
if (CtrlvsPCOS.flag == TRUE) {
  
  # Subset seurat object for only baseline
  Idents(object = x.seurat) = "Group_Stage"
  x.seurat = subset(x.seurat, idents = c("Control", "PCOS_W0"))
  
}

# Setting the idents and assay
Idents(object = x.seurat) <- paste0(celltype, "_labelled")
DefaultAssay(x.seurat) = "RNA"
Norm.assay = "RNA"

# Extract the celltypes to average
celltypes.seurat = levels(x.seurat)

for (celltype.x in celltypes.seurat) {
  
  # Generate a cell subtype output directory
  if (dir.exists(path = paste0("Output/8_DEG_GO_plotting_", celltype, "/", celltype.x)) == FALSE) {
    print(paste0("Generating output directory Output/8_DEG_GO_plotting_", celltype, "/", celltype.x))
    dir.create(path = paste0("Output/8_DEG_GO_plotting_", celltype, "/", celltype.x), recursive = TRUE)
    celltype.dir = paste0("Output/8_DEG_GO_plotting_", celltype, "/", celltype.x, "/")
  } else if (dir.exists(path = paste0("Output/8_DEG_GO_plotting_", celltype, "/", celltype.x)) == TRUE) {
    print("Directory exists")
    celltype.dir = paste0("Output/8_DEG_GO_plotting_", celltype, "/", celltype.x, "/")
  } else {
    print("Error with output directory")
  }
  
  # Subset seurat object for current celltype
  subset.cell = subset(x.seurat, idents = celltype.x)
  
  if (CtrlvsPCOS.flag == TRUE) {
    
    # Setting idents and assay on subset
    Idents(object = subset.cell) <- "orig.ident"
    DefaultAssay(subset.cell) = "RNA"
    
    # Do heatmap by calculating average expression
    subset.average = AverageExpression(subset.cell, return.seurat = TRUE, assays = "RNA")

    # Save the labelled and re-ordered object
    SaveH5Seurat(subset.cell, paste0(celltype.dir,Project_name, "_", celltype.x, "_average_CtrlvsPCOS.h5seurat"), overwrite = TRUE)
    
  } else if (CtrlvsPCOS.flag == FALSE) {
    
    # Setting idents and assay on subset
    Idents(object = subset.cell) <- "Group_Stage"
    DefaultAssay(subset.cell) = "RNA"
    
    # Do heatmap by calculating average expression
    subset.average = AverageExpression(subset.cell, return.seurat = TRUE, assays = "RNA")

    # Save the labelled and re-ordered object
    SaveH5Seurat(subset.average, paste0(celltype.dir,Project_name, "_", celltype.x, "_average_Group.h5seurat"), overwrite = TRUE)
    
  }
  
  
}