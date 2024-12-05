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
library(writexl)

# Setting the input files
Ref.group = "Labelled_Clusters_SCT.1" # "Labelled_Clusters_SCT.1"  #"Combined_labels"
Query.group = "seurat_clusters"
Transfer.type = "UMAP_integrated" # Transfer_data OR UMAP_integrated
Input.dir = "Output/1_Label_Transfer_Annotation/"
only.subsetting.flag = FALSE
Subset.percantage_proportion.threshold = 1 # Subset clusters with proportion higher than set. Default is 1%
set.CT_subset = "Immune" # Set to subset specific cell type


# Load the seurat objects
if (Transfer.type == "Transfer_data") {
  Query.list = list.files(Input.dir, pattern = "*Transfer_data_labelled.h5seurat", full.names = TRUE, recursive = TRUE)
  Query.list = grep(Ref.group, Query.list, value = TRUE)
  print(Query.list)
  # Predicted variable name
  seurat.ident = "predicted.id"
  seurat.predictQC = "prediction.score.max"
} else if (Transfer.type == "UMAP_integrated") {
  Query.list = list.files(Input.dir, pattern = "*UMAP_integrated_labelled.h5seurat", full.names = TRUE, recursive = TRUE)
  Query.list = grep(Ref.group, Query.list, value = TRUE)
  print(Query.list)
  # Predicted variable name
  seurat.ident = "predicted.celltype"
  seurat.predictQC = "predicted.celltype.score"
}

Automated_Annotation_plotting <- function(seurat.x, marker.list = Celltype.marker.list, 
                                          only.subsetting = only.subsetting.flag) {
  
  # Extract the query sample name
  seurat.sample = gsub(".*Output/1_Label_Transfer_Annotation//(.*?)/.*", "\\1", seurat.x)
  print(paste("Loading query seurat", seurat.sample))
  
  # Generate sample specific dircetorySetting the output directory
  if (dir.exists(path = paste0("Output/1_Label_Transfer_Annotation/", seurat.sample)) == FALSE) {
    print(paste0("Output/1_Label_Transfer_Annotation/", seurat.sample))
    dir.create(path =paste0("Output/1_Label_Transfer_Annotation/", seurat.sample), recursive = TRUE)
    Sample.dir = paste0("Output/1_Label_Transfer_Annotation/", seurat.sample, "/")
  } else if (dir.exists(path = "Output/1_Label_Transfer_Annotation") == TRUE) {
    print("Directory exists")
    Sample.dir = paste0("Output/1_Label_Transfer_Annotation/", seurat.sample, "/")
  } else {
    print("Error with output directory")
  }
  
  # Add transfer data type to Seurat sample name
  seurat.sample = paste0(seurat.sample, "_", Transfer.type)
  
  # Load the seurat object
  seurat.x = LoadH5Seurat(file = seurat.x)
  
  # Set idents
  Idents(seurat.x) = seurat.ident
  
  # Run FindAllMarkers
  DEG.res = FindAllMarkers(seurat.x, logfc.threshold = 0.1, min.pct = 0.1, assay = "RNA", )
  write_xlsx(DEG.res, path = paste0(Sample.dir, seurat.sample, "_FindAllMarkers_table.xlsx"))
  
  #Prepare the gene list
  marker.vector = unlist(marker.list)
  marker.vector = unique(marker.vector)
  
  # Filter the table
  DEG.res = DEG.res[DEG.res$gene %in% marker.vector,]
  
  # Add regulation column
  DEG.res$regulation = "NA"
  
  # Update the regulation column based on the avg_log2fc column
  DEG.res$regulation <- ifelse(DEG.res$avg_log2FC > 0, "positive", "negative")
  
  # Run FindAllMarkers with selected markers
  write_xlsx(DEG.res, path = paste0(Sample.dir, seurat.sample, "_FindAllMarkers_filtered_table.xlsx"))
  
  if (only.subsetting == FALSE) {
    
    # Generate UMAPs of the loaded query Seurat object
    DimPlot(seurat.x, reduction = "umap", group.by = seurat.ident)
    ggsave2(paste0(Sample.dir, seurat.sample, "_UMAP_predicted_celltype.pdf"))
    
    DimPlot(seurat.x, reduction = "umap", group.by = seurat.ident, label = TRUE)
    ggsave2(paste0(Sample.dir, seurat.sample, "_UMAP_labelled_predicted_celltype.pdf"))
    
    FeaturePlot(seurat.x, features = seurat.predictQC, reduction = "umap", pt.size = 0.1)
    ggsave2(paste0(Sample.dir, seurat.sample, "_UMAP_predicted_celltype_score.pdf"))
    
    DimPlot(seurat.x, reduction = "spatial", group.by = seurat.ident)
    ggsave2(paste0(Sample.dir, seurat.sample, "_Spatial_predicted_celltype.pdf"))
    
    FeaturePlot(seurat.x, features = seurat.predictQC, reduction = "spatial")
    ggsave2(paste0(Sample.dir, seurat.sample, "_Spatial_predicted_celltype_score.pdf"))
    
    FeaturePlot(seurat.x, features = seurat.predictQC, reduction = "umap")
    ggsave2(paste0(Sample.dir, seurat.sample, "_UMAP_predicted_celltype_score.pdf"))
    
    # Check marker expression
    for (marker.x in names(marker.list)) {
      
      print(marker.x)
      marker.plot = marker.list[[marker.x]]
      
      DotPlot(object = seurat.x, features = marker.plot, cluster.idents = TRUE, scale = FALSE) + RotatedAxis()
      ggsave2(paste0(Sample.dir, marker.x, "_", seurat.sample, "_Marker_Dotplot.pdf"))
      
    }
    
  }
  
  # Check and plot proportions of celltypes
  prop.df = data.frame(prop.table(table(seurat.x[[seurat.ident]]))*100)
  colnames(prop.df) = c("Celltype", "Percentage")
  write_xlsx(prop.df, path = paste0(Sample.dir, seurat.sample, "_Celltype_proportion_table.xlsx"))
  
  ggplot(prop.df, aes(x = "", y = Percentage, fill = Celltype)) +
    geom_bar(stat = "identity") + theme_cowplot() +
    theme(legend.position = "right") +
    labs(x = "", y = "Percentage", fill = "Cell Type")
  
  ggsave2(paste0(Sample.dir, seurat.sample, "_barplot_celltypes.pdf"), 
          dpi = 700)
  
  # Generate sample specific dircetorySetting the output directory of subset celltypes
  if (dir.exists(path = paste0("Output/2_Automated_Subset_Seurat/", seurat.sample)) == FALSE) {
    print(paste0("Output/2_Automated_Subset_Seurat/", seurat.sample))
    dir.create(path =paste0("Output/2_Automated_Subset_Seurat/", seurat.sample), recursive = TRUE)
    Subset.dir = paste0("Output/2_Automated_Subset_Seurat/", seurat.sample, "/")
  } else if (dir.exists(path = "Output/2_Automated_Subset_Seurat") == TRUE) {
    print("Directory exists")
    Subset.dir = paste0("Output/2_Automated_Subset_Seurat/", seurat.sample, "/")
  } else {
    print("Error with output directory")
  }
  
  # Subset the annotated celltypes. Annotate only celltype with proportions higher than 1%.
  # Small subtypes generate errors
  
  if (is.null(set.CT_subset) == TRUE) {
    CT.subsets = prop.df$Celltype[prop.df$Percentage > Subset.percantage_proportion.threshold]
    print(CT.subsets)
    print(paste0("Subsetting celltype with proportion higher than ", Subset.percantage_proportion.threshold, "% ", CT.subsets))
  } else if (is.null(set.CT_subset) == FALSE) {
    CT.subsets = set.CT_subset
    print(paste0("Subsetting selected celltype ", CT.subsets))
    
  }
  
  for (CT.x in CT.subsets) {
    
    # Check if subsetting has already been done
    CT.check = list.files(Subset.dir, pattern = paste0(".*", Transfer.type, "_", CT.x, "_subset.h5Seurat"), full.names = FALSE)
    
    if (length(CT.check) > 0) {
      print(paste("Subsetting has already been performed on", CT.check))
      next()
    } else if (length(CT.check) == 0) {
      # Subset the seurat object
      print(paste("Subsetting", CT.x))
      CT.subset = subset(seurat.x, idents = CT.x)
    }
    
    # Subset the seurat object
    print(paste("Subsetting", CT.x))
    CT.subset = subset(seurat.x, idents = CT.x)
    
    # Remove seurat object if last CT is subseted. Do this to reduce memory.
    if (CT.subsets[length(CT.subsets)] == CT.x) {
      print("Last sample to subset, removing seurat object")
      seurat.x = NULL
    }
    
    # Run the subset analysis function
    Do.subsetting = CT.subsetting(CT.subset, CT.sample = CT.x, 
                                  output.subset = paste0(Subset.dir, seurat.sample))
    
  }
  
  # Run the subsetting function
  #Do.subsetting = lapply(CT.subsets, CT.subsetting, output.subset = paste0(Subset.dir, seurat.sample))
  
}

# Subsetting function
CT.subsetting <- function(CT.subset, CT.sample = CT.x, output.subset) {
  
  # Normalise the subset object and scale
  DefaultAssay(CT.subset) <- "RNA"
  CT.subset = NormalizeData(CT.subset, normalization.method = "LogNormalize")
  CT.subset = FindVariableFeatures(CT.subset, selection.method = "vst", verbose = T)
  top10_var <- head(VariableFeatures(CT.subset), 10)
  all.genes = rownames(CT.subset)
  CT.subset = ScaleData(CT.subset, verbose = T, features = all.genes)
  
  # Run scTransform and cluster the data
  # Perform scTransformation normalisation
  DefaultAssay(CT.subset) <- "RNA"
  CT.subset = SCTransform(CT.subset, verbose = T)
  
  # Do PCA on scTranformed object
  CT.subset = RunPCA(CT.subset)
  
  # Plotting PCA results
  PCA_plot = DimPlot(CT.subset, reduction = "pca")
  ggsave2(plot = PCA_plot, filename = paste0(output.subset, "_", CT.sample, "_scTransform_PCA_plot.pdf"))
  
  PCA_heatmap = DimHeatmap(CT.subset, dims = 1:15, cells = 500, balanced = TRUE)
  ggsave2(plot = PCA_heatmap, filename = paste0(output.subset, "_", CT.sample, "_scTransform_PCA_heatmap.pdf"))
  
  # Determine he dimensionality of the data
  elbow_plot = ElbowPlot(CT.subset, reduction = "pca", ndims = 25)
  ggsave2(plot = elbow_plot, filename = paste0(output.subset, "_", CT.sample, "_scTransform_elbow_plot.pdf"))
  
  # Do UMAP on seurat object
  CT.subset = RunUMAP(CT.subset, dims = 1:30)
  
  # Run clustering on the object
  CT.subset = FindNeighbors(object = CT.subset, reduction = "pca", dims = 1:30, verbose = FALSE)
  CT.subset = FindClusters(object = CT.subset, resolution = 0.7, verbose = FALSE)
  
  # Plotting UMAP
  UMAP_plot = DimPlot(CT.subset, label = TRUE, repel = TRUE)
  ggsave2(plot = UMAP_plot, filename = paste0(output.subset, "_", CT.sample, "_scTransform_UMAP.pdf"))
  
  # Plotting UMAP on spatial
  UMAP_plot = DimPlot(CT.subset, label = FALSE, repel = FALSE, reduction = "spatial")
  ggsave2(plot = UMAP_plot, filename = paste0(output.subset, "_", CT.sample, "_scTransform_UMAP_spatial.pdf"))
  
  # Reset deafault to RNA
  DefaultAssay(CT.subset) <- "RNA"
  
  # Saved the subset Seurat object
  SaveH5Seurat(object = CT.subset, filename = paste0(output.subset, "_", CT.sample, "_subset.h5Seurat"), overwrite = TRUE)
  
  # Set subset to NULL to save memory
  CT.subset = NULL
  
  print(paste("Finished subsetting", CT.sample))
  
}

# Marker list
Celltype.marker.list = list()
Celltype.marker.list[["Major_markers"]] = c("EPCAM", "CPM", "LGR5", "IGF1", "DCN", "COL6A1", "GUCY1A2", "ACTA2", "NOTCH3", "CD14",
                                            "CSF1R", "LYZ", "STK17B", "NCAM1", "CCL5", "CD2", "PCDH17", "VWF", "PROX1", "FLT4")
Celltype.marker.list[["Epithelium_Subtypes"]] = c('PTGS1', 'VTCN1', 'SLC26A7', 'LGR5', 'KRT5', 'WNT7A', 'CPM', 
                                                    'IHH', 'EMID1', 'PPARG', 'C2CD4A', 'SLC18A2','PAEP', 'CXCL14', 
                                                    'MKI67', 'HMGB2', 'AR', 'CDC20B', 'CCNO', 'HES6')
Celltype.marker.list[["Stroma_Subtypes"]] = c('ESR1', 'PGR', 'IGF1', 'ECM1', 'PAEP', 'OGN', 'TOP2A', 'MKI67', 
                                            'THY1', 'COL1A1', 'PCOLCE', 'C7', 'ACTA2', 'ACTG2', 'MCAM')
Celltype.marker.list[["Endothelial_Lymphatic_Subtypes"]] = c("PECAM1", "CD34" ,
                                                           "ACKR1", "PLVAP", # Endothelial Vein
                                                           "SEMA3G", "GJA5", # Endothelial Artery
                                                           "TOP2A", "MKI67", # Endothelial proliferative
                                                           "COL3A1", "WNT5A", "MMP11", # Mesenchymal
                                                           "PROX1", "FLT4") # Lymphatic
Celltype.marker.list[["Immune_Subtypes"]] = c("FOXP3", "CD3G", "CD27", "IKZF4", "CD8A", "CCL5", "CD4", 
                                            "IL7R", "RORC", #Treg, T-cell CD8, T-cell CD4, ILC3
                                            "MS4A1", "IGHM", #B-celler
                                            "NCAM1", "ITGA1", "SPINK2", "CSF1", "CD160", "GNLY", # uNK 1-3
                                            "IL3RA", "LILRA4", "PLD4",# pDC
                                            "EBI3", "CCR7", "CCL19", # Migratory DC
                                            "BATF3", "CADM1", "CLEC9A", #DC1
                                            "CLEC10A", "FCER1A", "CD1C", #DC2
                                            "CD14", "SELENOP", "HMOX1", "IL1B", #uM 1-2
                                            "CPA3", "KIT", "MS4A2") # Mast cells

# Run the plotting function
plot.res = lapply(Query.list, Automated_Annotation_plotting)
