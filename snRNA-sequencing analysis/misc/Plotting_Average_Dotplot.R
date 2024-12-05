#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(pheatmap)

## Set variables
Project.name = "Endo_Control_PCOS"
Group.cols = c("#A0A0A0", "#D098B2", "#95BFE1", "#65D46E")
Group.order = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS")
target.features = c("ESR1", "PGR", "AR")
CT.order = c("Lumenal", "SOX9+ LGR5+", "SOX9+ LGR5-", "SOX9+ proliferative", "AR+", "Ciliated", # Epithelium
             "Stroma 1", "Stroma 2", "Stroma proliferative", "Fibroblast", "uSMC", # Stromal
             "uM 1", "uM 2", "T-cells CD4+", "T-cells CD8+", "uNK 1", "uNK 2", "uNK 3", # Immune
             "Endothelial Artery", "Endothelial Vein", "Lymphatic") # Endothelium and lymphatic


# Set the input Seurat
seurat.path = "Data/Seurat_Subset/Endo_Combined_Control_PCOS_W0_subset.h5seurat"

# Generate or set the output directory
if (dir.exists(path = paste0("Output/Plotting_Average_Dotplot")) == FALSE) {
  print(paste0("Generating output directory Output/Plotting_Average_Dotplot"))
  dir.create(path = paste0("Output/Plotting_Average_Dotplot"), recursive = TRUE)
  Output.dir = paste0("Output/Plotting_Average_Dotplot/")
} else if (dir.exists(path = paste0("Output/Plotting_Average_Dotplot")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/Plotting_Average_Dotplot/")
} else {
  print("Error with output directory")
}

# Set the output name
output.path = paste0(Output.dir, Project.name)

# Load the seurat object, only RNA to decrease memory load
seurat.x = LoadH5Seurat(file = seurat.path, assay = "RNA")

# Subset selected celltypes
Idents(object = seurat.x) <- "Combined_labels"
print(levels(seurat.x))
seurat.x = subset(x = seurat.x, idents = CT.order)
#Idents(seurat.x) <- factor(Idents(seurat.x), levels = rev(CT.order))
seurat.x$Combined_labels = factor(seurat.x$Combined_labels, levels = rev(CT.order))
#seurat.x$Combined_labels = factor(seurat.x$Combined_labels, levels = rev(CT.order))
print(levels(seurat.x$Combined_labels))

# Normalize the data
seurat.x = NormalizeData(seurat.x)

# Setting the order of the dotplot
# Step 1: Save the original levels of the cell types
levels.order <- levels(seurat.x$Combined_labels)
# Step 2: Create the combined group by combining the original levels with the condition
seurat.x$Combined_labels_Group <- paste(seurat.x$Combined_labels, seurat.x$Group_Stage, sep = "_")
# Step 3: Create a new set of levels by interleaving the groups based on cell type
group.levels <- as.vector(sapply(levels.order, function(x) c(paste(x, "PCOS_W0", sep = "_"), paste(x, "Control", sep = "_"))))
# Step 4: Set the levels for Combined_labels_Group to enforce the desired order
seurat.x$Combined_labels_Group <- factor(seurat.x$Combined_labels_Group, levels = group.levels)
# Check the new order
print(levels(seurat.x$Combined_labels_Group))

Idents(object = seurat.x) <- "Combined_labels_Group"
print(levels(seurat.x))

plot.x = DotPlot(object = seurat.x, assay = "RNA", features = target.features, scale = FALSE, dot.scale = 9,
        cols = c("lightgrey", "#CC0000"))
ggsave(filename = paste0(output.path, "_Dotplot_Combined_Labels.pdf"),
       plot = plot.x, dpi = 300, width = 7, height = 15)


DotPlot(object = seurat.x, features = target.features, assay = "RNA", cols =c(), scale = TRUE, dot.scale = 8,
        group.by = "Group_Stage", split.by = "Combined_labels_Gro") +
  scale_colour_gradient() +
  guides(color = guide_colorbar(title = 'Average Expression'))
