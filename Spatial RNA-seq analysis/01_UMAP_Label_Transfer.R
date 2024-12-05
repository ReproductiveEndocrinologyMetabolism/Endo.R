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

# Set parameters
Ref.name = "Endo_10x"
Input.dir = "Output/0_Preprocessing/"
load.anchors = FALSE
integration.feat = 2500
Transfer.data.flag = FALSE
UMAP.integration.flag = TRUE
Ref.seurat = "/proj/nobackup/sens2022003/EndoR_220925_10%mtDNA_filter-mtgenes_no-ref_2500feat/Output/9_Main_relabelled/Endo_All_Combined_labels.h5seurat"
Ref.group = c("Combined_labels") # "Labelled_Clusters_SCT.1"  #"Combined_labels"
Query.group = "seurat_clusters"

# Load the seurat objects
Query.list = list.files(Input.dir, pattern = "*.h5Seurat", full.names = TRUE, recursive = TRUE)
print(Query.list)

# Setting the output directory
if (dir.exists(path = "Output/1_Full_Integration_Automated_Annotation") == FALSE) {
  print("Output/1_Full_Integration_Automated_Annotation")
  dir.create(path = "Output/1_Full_Integration_Automated_Annotation", recursive = TRUE)
  Output.dir = "Output/1_Full_Integration_Automated_Annotation/"
} else if (dir.exists(path = "Output/1_Full_Integration_Automated_Annotation") == TRUE) {
  print("Directory exists")
  Output.dir = "Output/1_Full_Integration_Automated_Annotation/"
} else {
  print("Error with output directory")
}

# Run the integration of the individually filtered spatial object
seurat.list = lapply(Query.list, function(x) LoadH5Seurat(file = x))

print("Integration of SCT transformed object")
# Preparing SC transformed objects for integration
seurat.list <- lapply(X = seurat.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = integration.feat)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
seurat.list <- lapply(X = seurat.list, FUN = RunPCA, features = features)

print("Finding integration anchors")
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                         normalization.method = "SCT", 
                                         anchor.features = features,
                                         dims = 1:30,
                                         reduction = "rpca", 
                                         k.anchor = 5)

# Remove the seurat.list to reduce memory load
seurat.list = NULL

# The anchor dataset is saved
saveRDS(seurat.anchors, paste0(Output.dir, "Spatial_integration_seurat_anchors.RDS"))

# Data integration
print("Integrating the data")
seurat.integrated <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", dims = 1:30)

# Remove seurat.anchors
seurat.anchors = NULL

# Cluster the integrated data
seurat.integrated <- RunPCA(seurat.integrated, verbose = TRUE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30)
seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30)
seurat.integrated <- FindClusters(seurat.integrated)

# Saving the seurat object
print("Saving the integrated data")
SaveH5Seurat(object = seurat.integrated, filename = paste0(Output.dir,"Spatial_filtered_integrated.h5seurat"), overwrite = TRUE)

# Log2 normalising 
print("Normalising the integrated data")
DefaultAssay(seurat.integrated) = "RNA"
seurat.integrated = NormalizeData(seurat.integrated)
seurat.integrated = FindVariableFeatures(seurat.integrated)
all.genes = rownames(seurat.integrated)
seurat.integrated = ScaleData(seurat.integrated, features = all.genes)

# Saving the seurat object
SaveH5Seurat(object = seurat.integrated, filename = paste0(Output.dir,"Spatial_filtered_integrated.h5seurat"), overwrite = TRUE)

# Renaming for groups and stage
print("Plotting the integrated data")
DimPlot(seurat.integrated, reduction = "umap", group.by = "orig.ident")
ggsave2(paste0(Output.dir, "Spatial_filtered_integrated_orig_ident_UMAP.pdf"), dpi = 300)

# Renaming for groups and stage
DimPlot(seurat.integrated, reduction = "umap", group.by = "seurat_clusters")
ggsave2(paste0(Output.dir, "Spatial_filtered_integrated_seurat_clusters_UMAP.pdf"), dpi = 300)

# Load the Reference and Query seurat
print("Loading reference seurat")
Ref.seurat = LoadH5Seurat(file = Ref.seurat)

# Generate UMAPs of the loaded reference Seurat object
DimPlot(Ref.seurat, reduction = "umap", group.by = Ref.group)
ggsave2(paste0(Output.dir, "Reference_", Ref.group, "_UMAP.pdf"), dpi = 300)

# Setting assays to RNA for log2 normalisation
DefaultAssay(Ref.seurat) = "RNA"

# Check assay names so that they match
print("Checking assay for referemce:")
Assays(Ref.seurat)

# Set idents
Idents(Ref.seurat) = Ref.group

# Perform standard preprocessing on reference to detect variable features
Ref.seurat = NormalizeData(Ref.seurat)
Ref.seurat = FindVariableFeatures(Ref.seurat)
Ref.seurat = ScaleData(Ref.seurat)

# Check common genes
common.genes = intersect(rownames(Ref.seurat), rownames(seurat.integrated))
head(common.genes)
print(paste("Number of common genes are:", length(common.genes)))

# Find Transfer Anchors between the seurat object
Anno.anchors = FindTransferAnchors(reference = Ref.seurat, query = seurat.integrated, dims = 1:30, reference.reduction = "pca")

# Do unimodal UMAP projection THIS HAS BEEN DONE BEFORE?
Ref.seurat <- RunUMAP(Ref.seurat, dims = 1:30, return.model = TRUE)

# Run MapQuery wrapper
seurat.integrated <- MapQuery(anchorset = Anno.anchors, reference = Ref.seurat, query = seurat.integrated,
                         refdata = list(celltype = Ref.group), reference.reduction = "pca", reduction.model = "umap")

print("Annotation finished, saving the labelled integrated seurat")
# Saving the seurat object
SaveH5Seurat(object = seurat.integrated, filename = paste0(Output.dir,"Spatial_filtered_integrated.h5seurat"), overwrite = TRUE)

# Remove the ref.seurat and anno.anchors to reduce memory load
Ref.seurat = NULL
Anno.anchors = NULL

print("Extracting annotation from each sample of the integrated object")

# Extracting predicted celltypes of each sample from the integrated object
# Set the ident with selected label
Idents(seurat.integrated) = "orig.ident"

for (Query.sample in Query.list) {
  
  # Extract the query sample name
  seurat.sample = gsub(".*Output/0_Preprocessing//(.*?)/.*", "\\1", Query.sample)
  
  # Subset the integrated seurat
  seurat.subset = subset(seurat.integrated, subset = orig.ident == seurat.sample)
  
  # Extract the predicted labels
  seurat.label = seurat.subset$predicted.celltype
  colnames(seurat.label) = "Predicted_integrated_label"
  seurat.subset = NULL
  
  # Save the celltype labels
  saveRDS(seurat.label, file = paste0(Output.dir, seurat.sample, "_Labels.rds"))

  # Load the query sample to project the labels on
  seurat.x = LoadH5Seurat(Query.sample)
  
  seurat.x$Predicted_Integrated_Labels = seurat.label
  
  # Save the object
  SaveH5Seurat(object = seurat.x, filename = paste0(Output.dir, seurat.sample,  "_Integrated_Annotation_filtered.h5seurat"), overwrite = TRUE)
  
  # Generate UMAPs of the loaded reference Seurat object
  DimPlot(seurat.x, reduction = "umap", group.by = "Predicted_Integrated_Labels")
  ggsave2(paste0(Output.dir, seurat.sample, "_Predicted_Integrated_Labels_UMAP"), dpi = 300)
  
  # Generate UMAPs of the loaded reference Seurat object
  DimPlot(seurat.x, reduction = "spatial", group.by = "Predicted_Integrated_Labels")
  ggsave2(paste0(Output.dir, seurat.sample, "_Predicted_Integrated_Labels_Spatial.pdf"), dpi = 300)
  
}


