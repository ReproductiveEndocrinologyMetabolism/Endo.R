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
library(readxl)
library(writexl)

# Set parameters
seurat.ident = "predicted.celltype" # "seurat_cluster" or "predicted.celltype"
Norm.assay = "RNA"
Transfer.type = "UMAP_integrated" # UMAP_integrated, Transfer_data or "" for manual
Input.dir = "Output/2_Automated_Subset_Seurat/"  #"Output/2_Automated_Subset_Seurat/" or "Output/2_Manual_Subset_Seurat/" 
                                                  #or "Output/1_Label_Transfer_Annotation/"
annotations.set = FALSE

# Load the seurat objects
Epithelium.list = list.files(Input.dir, pattern = "*Epithelium_subset.h5Seurat", full.names = TRUE, recursive = TRUE)
Epithelium.list = grep(Transfer.type, Epithelium.list, value = TRUE)
Stroma.list = list.files(Input.dir, pattern = "*Stroma_subset.h5Seurat", full.names = TRUE, recursive = TRUE)
Stroma.list = grep(Transfer.type, Stroma.list, value = TRUE)
Immune.list = list.files(Input.dir, pattern = "*Immune_subset.h5Seurat", full.names = TRUE, recursive = TRUE)
Immune.list = grep(Transfer.type, Immune.list, value = TRUE)
Query.group = "seurat_clusters"

# Load excel with annotations after setting these
if (annotations.set == TRUE) {
  epithelium.table = as.data.frame(read_excel("Data/Epithelium_cluster_annotation.xlsx"))
  stroma.table = as.data.frame(read_excel("Data/Stroma_cluster_annotation.xlsx"))
  immune.table = as.data.frame(read_excel("Data/Immune_cluster_annotation.xlsx"))
  endothelial.table = as.data.frame(read_excel("Data/Endothelial_cluster_annotation.xlsx"))
}

Subset_marker_plotting <- function(seurat.x, marker.list = Celltype.marker.list, 
                                   cluster.table = NULL, CT.sample = NULL) {
  
  # Extract subset type
  #CT.sample = gsub(paste0(".*_bin30_", Transfer.type, "_([A-Za-z]+)_subset\\.h5Seurat"), "\\1", seurat.x)
  
  # Extract the query sample name
  seurat.sample = gsub(".*Output/2_.*_Subset_Seurat//(.*?)/.*", "\\1", seurat.x)
  
  # Generate output directory
  if (dir.exists(path = paste0("Output/3_Subtype_Annotation/", seurat.sample)) == FALSE) {
    print(paste0("Output/3_Subtype_Annotation/", seurat.sample))
    dir.create(path =paste0("Output/3_Subtype_Annotation/", seurat.sample), recursive = TRUE)
    Output.dir = paste0("Output/3_Subtype_Annotation/", seurat.sample, "/")
  } else if (dir.exists(path = paste0("Output/3_Subtype_Annotation/", seurat.sample)) == TRUE) {
    print("Directory exists")
    Output.dir = paste0("Output/3_Subtype_Annotation/", seurat.sample, "/")
  } else {
    print("Error with output directory")
  }
  
  # Load the seurat object
  seurat.x = LoadH5Seurat(file = seurat.x)
  
  # Set idents
  Idents(seurat.x) = seurat.ident
  
  if (annotations.set == FALSE) {
    
    # Add subset information
    #seurat.sample = paste0(seurat.sample, "_", CT.sample, "_subset")
    print(paste("Plotting subset", seurat.sample))
    
    DimPlot(seurat.x, reduction = "umap", label = TRUE)
    ggsave2(paste0(Output.dir, CT.sample, "_Seurat_cluster_UMAP.pdf"))
    
    DimPlot(seurat.x, reduction = "spatial", label = TRUE)
    ggsave2(paste0(Output.dir, CT.sample, "_Seurat_cluster_spatial.pdf"))
    
    # Plot the markers
    marker.x = names(marker.list)[1]
    for (marker.x in names(marker.list)) {
      
      print(marker.x)
      marker.plot = marker.list[[marker.x]]
      
      if (length(levels(seurat.x[[seurat.ident]])) > 1) {
        DotPlot(object = seurat.x, features = marker.plot, cluster.idents = TRUE, scale = FALSE) + RotatedAxis()
        ggsave2(paste0(Output.dir, CT.sample, "_", marker.x, "_Marker_Dotplot.pdf"))
      }
      
      VlnPlot(seurat.x, features = marker.plot, sort = TRUE, pt.size = 0.5)
      ggsave2(paste0(Output.dir, CT.sample, "_", marker.x, "_Marker_Vlnplot.pdf"))
      
      VlnPlot(seurat.x, features = marker.plot, sort = TRUE, pt.size = 0.5)
      ggsave2(paste0(Output.dir, CT.sample, "_", marker.x, "_Marker_Vlnplot_small.pdf"), height = 3, width = 3)
      
      FeaturePlot(seurat.x, features = marker.plot, reduction = "umap", pt.size = 0.5, label = FALSE, label.size = 2)
      ggsave2(paste0(Output.dir, CT.sample, "_", marker.x, "_Marker_Featureplot_UMAP.pdf"), height = 3, width = 6)
      
                  
    }
    
    # Plot additional features
    # Plot the cell cycle score
    # Visualize the distribution of cell cycle markers across
    VlnPlot(seurat.x, features = c("S.Score","G2M.Score"), group.by = seurat.ident, pt.size = 0)
    ggsave2(paste0(Output.dir, CT.sample, "_Cell.Cycle_Phase_VlnPlot.pdf"))
    
    # Plot cell cycle phases
    Phase_UMAP = DimPlot(seurat.x, group.by = "Phase", reduction = "umap")
    ggsave2(plot = Phase_UMAP, filename = paste0(Output.dir, CT.sample, "_Cell-Cycle_Phase_UMAP_filtered.pdf"))
    
    Phase_Spatial = DimPlot(seurat.x, group.by = "Phase", reduction = "spatial")
    ggsave2(plot = Phase_Spatial, filename = paste0(Output.dir, CT.sample, "_Cell-Cycle_Phase_spatial_filtered.pdf"))
    
    # Plotting additional UMAPs before saving
    FeaturePlot(seurat.x, features = "percent.mt", reduction = "umap", label = TRUE)
    ggsave2(paste0(Output.dir, CT.sample, "_UMAP_percent_mt.pdf"))
    
    VlnPlot(seurat.x, features = "percent.mt", sort = TRUE, pt.size = 0)
    ggsave2(paste0(Output.dir, CT.sample, "_Vlnplot_percent_mt.pdf"))
    
    FeaturePlot(seurat.x, features = "percent.hb", reduction = "umap", label = TRUE)
    ggsave2(paste0(Output.dir, CT.sample, "_UMAP_percent_hb.pdf"))
    
    FeaturePlot(seurat.x, features = "S.Score", reduction = "umap")
    ggsave2(paste0(Output.dir, CT.sample, "_UMAP_S_score.pdf"))
    
    FeaturePlot(seurat.x, features = "G2M.Score", reduction = "umap")
    ggsave2(paste0(Output.dir, CT.sample, "_UMAP_G2M_Score.pdf"))
    
  } else if (annotations.set == TRUE) {
    
    # Add subset information
    seurat.sample = paste0(seurat.sample, "_", CT.sample, "_subset")
    print(paste("Loading query seurat", seurat.sample))
    
    # Extracting cluster label information
    cluster.table = cluster.table[,c("Cluster", seurat.sample)]
    cluster.table = na.omit(cluster.table)
    cluster.table = cluster.table[,seurat.sample]
    
    # Set the labels in the Seurat object
    names(cluster.table) = levels(seurat.x)
    seurat.x <- RenameIdents(seurat.x, cluster.table)
    seurat.x[["Labelled_subcluster"]] <- Idents(seurat.x)
    
    # Plot the labelled object
    plot.dim = DimPlot(seurat.x, reduction = "spatial", label = FALSE, pt.size = 0.5)
    ggsave2(paste0(Output.dir, CT.sample, "_Labelled_Clusters_Spatial.pdf"))
    plot.dim = DimPlot(seurat.x, reduction = "umap", label = FALSE, pt.size = 0.5)
    ggsave2(paste0(Output.dir, CT.sample, "_Labelled_Clusters_umap.pdf"))
    
    # Save the labelled subset
    SaveH5Seurat(object = seurat.x, filename = paste0(Output.dir, seurat.sample, "_labelled.h5Seurat"), overwrite = TRUE)
    
    
    # Plot the markers on the annotated subcluster
    for (marker.x in names(marker.list)) {
      
      print(marker.x)
      marker.plot = marker.list[[marker.x]]
      
      DotPlot(object = seurat.x, features = marker.plot, cluster.idents = TRUE, scale = FALSE) + RotatedAxis()
      ggsave2(paste0(Output.dir, CT.sample, "_", marker.x,"_Marker_Dotplot_labelled.pdf"))
      
      VlnPlot(seurat.x, features = marker.plot, sort = TRUE, pt.size = 0.01)
      ggsave2(paste0(Output.dir, CT.sample, "_", marker.x, "_Marker_Vlnplot_labelled.pdf"))
      
      FeaturePlot(seurat.x, features = marker.plot, pt.size = 0.001, label = TRUE, label.size = 2)
      ggsave2(paste0(Output.dir, marker.x, "_Marker_Featureplot_UMAP_labelled.pdf"))
      
    }
    
  }
  
}

# Marker list
# Updated cell cluster markers
Celltype.marker.list = list()
Celltype.marker.list[["Stroma_Subtypes"]] = c('ESR1', 'PGR', 'IGF1', 'ECM1', 'PAEP', 'OGN', 'TOP2A', 'MKI67', 
                                              'THY1', 'COL1A1', 'PCOLCE', 'C7', 'ACTA2', 'ACTG2', 'MCAM')
Celltype.marker.list[["Epithelium_Subtypes"]] = c('PTGS1', 'VTCN1', 'SLC26A7', 'LGR5', 'KRT5', 'WNT7A', 'CPM', 
                                                  'IHH', 'EMID1', 'PPARG', 'C2CD4A', 'SLC18A2','PAEP', 'CXCL14', 
                                                  'MKI67', 'HMGB2', 'AR', 'CDC20B', 'CCNO', 'HES6')
Celltype.marker.list[["Immune_Subtypes"]] = c("FOXP3", "CD3G", "CD27", "IKZF4", "CD8A", "CCL5", "CD4", 
                                              "IL7R", "RORC", #Treg, T-cell CD8, T-cell CD4, ILC3
                                              "MS4A1", "IGHM", #B-celler
                                              "NCAM1", "ITGA1", "SPINK2", "CSF1", "CD160", "GNLY", # uNK 1-3
                                              "IL3RA", "LILRA4", "PLD4",# pDC
                                              "EBI3", "CCR7", "CCL19", # Migratory DC
                                              "BATF3", "CADM1", "CLEC9A", #DC1
                                              "CLEC10A", "FCER1A", "CD1C", #DC2
                                              "CD14", "SELENOP", "HMOX1", "IL1B", #uM 1-2
                                              "CPA3", "KIT", "MS4A2", # Mast cells
                                              "PAEP") # Gene of interest
Celltype.marker.list[["Major_markers"]] = c("EPCAM", "CPM", "LGR5", "IGF1", "DCN", "COL6A1", "GUCY1A2", "ACTA2", "NOTCH3", "CD14",
                                            "CSF1R", "LYZ", "STK17B", "NCAM1", "CCL5", "CD2", "PCDH17", "VWF", "PROX1", "FLT4")

Epithelium.marker.list = list()
Epithelium.marker.list[["Epithelium_Subtypes"]] = c('PTGS1', 'VTCN1', 'SLC26A7', 'LGR5', 'KRT5', 'WNT7A', 'CPM', 
                                                    'IHH', 'EMID1', 'PPARG', 'C2CD4A', 'SLC18A2','PAEP', 'CXCL14', 
                                                    'MKI67', 'HMGB2', 'AR', 'CDC20B', 'CCNO', 'HES6')
Epithelium.marker.list[["Major_markers"]] = c("EPCAM", "CPM", "LGR5", "IGF1", "DCN", "COL6A1", "GUCY1A2", "ACTA2", "NOTCH3", "CD14",
                                              "CSF1R", "LYZ", "STK17B", "NCAM1", "CCL5", "CD2", "PCDH17", "VWF", "PROX1", "FLT4")

Stroma.marker.list = list()
Stroma.marker.list[["Stroma_Subtypes"]] = c('ESR1', 'PGR', 'IGF1', 'ECM1', 'PAEP', 'OGN', 'TOP2A', 'MKI67', 
                                            'THY1', 'COL1A1', 'PCOLCE', 'C7', 'ACTA2', 'ACTG2', 'MCAM')
Stroma.marker.list[["Endothelial_Lymphatic_Subtypes"]] = c("PECAM1", "CD34" ,
                                                           "ACKR1", "PLVAP", # Endothelial Vein
                                                           "SEMA3G", "GJA5", # Endothelial Artery
                                                           "TOP2A", "MKI67", # Endothelial proliferative
                                                           "COL3A1", "WNT5A", "MMP11", # Mesenchymal
                                                           "PROX1", "FLT4") # Lymphatic
Stroma.marker.list[["Immune_Subtypes"]] = c("FOXP3", "CD3G", "CD27", "IKZF4", "CD8A", "CCL5", "CD4", 
                                            "IL7R", "RORC", #Treg, T-cell CD8, T-cell CD4, ILC3
                                            "MS4A1", "IGHM", #B-celler
                                            "NCAM1", "ITGA1", "SPINK2", "CSF1", "CD160", "GNLY", # uNK 1-3
                                            "IL3RA", "LILRA4", "PLD4",# pDC
                                            "EBI3", "CCR7", "CCL19", # Migratory DC
                                            "BATF3", "CADM1", "CLEC9A", #DC1
                                            "CLEC10A", "FCER1A", "CD1C", #DC2
                                            "CD14", "SELENOP", "HMOX1", "IL1B", #uM 1-2
                                            "CPA3", "KIT", "MS4A2") # Mast cells
Stroma.marker.list[["Major_markers"]] = c("EPCAM", "CPM", "LGR5", "IGF1", "DCN", "COL6A1", "GUCY1A2", "ACTA2", "NOTCH3", "CD14",
                                          "CSF1R", "LYZ", "STK17B", "NCAM1", "CCL5", "CD2", "PCDH17", "VWF", "PROX1", "FLT4")

Immune.marker.list = list()
Immune.marker.list[["Immune_Subtypes"]] = c("FOXP3", "CD3G", "CD27", "IKZF4", "CD8A", "CCL5", "CD4", 
                                             "IL7R", "RORC", #Treg, T-cell CD8, T-cell CD4, ILC3
                                             "MS4A1", "IGHM", #B-celler
                                             "NCAM1", "ITGA1", "SPINK2", "CSF1", "CD160", "GNLY", # uNK 1-3
                                             "IL3RA", "LILRA4", "PLD4",# pDC
                                             "EBI3", "CCR7", "CCL19", # Migratory DC
                                             "BATF3", "CADM1", "CLEC9A", #DC1
                                             "CLEC10A", "FCER1A", "CD1C", #DC2
                                             "CD14", "SELENOP", "HMOX1", "IL1B", #uM 1-2
                                             "CPA3", "KIT", "MS4A2") # Mast cells
Immune.marker.list[["PAEP_uM_uNK_expression"]] = c("NCAM1", "ITGA1", "SPINK2", "CSF1", "CD160", "GNLY", # uNK 1-3
                                                   "CD14", "SELENOP", "HMOX1", "IL1B", #uM 1-2
                                                   "PAEP") #PAEP
Immune.marker.list[["PAEP_Immune_expression"]] = c("CD14","CSF1R", "LYZ", #Myeloid
                                                   "STK17B", "NCAM1", "CCL5", "CD2", #Lymphoid
                                                   "PAEP") #PAEP
Immune.marker.list[["PAEP_PTPRC_expression"]] = c("PTPRC", #Immune marker
                                                   "PAEP") #PAEP

# Running marker plotting
print("Plotting marker genes on subtypes")
plotting.res <- lapply(Epithelium.list, function(seurat.x) {
  Subset_marker_plotting(seurat.x, marker.list = Epithelium.marker.list, CT.sample = "Epithelium")
})
plotting.res <- lapply(Stroma.list, function(seurat.x) {
  Subset_marker_plotting(seurat.x, marker.list = Stroma.marker.list, CT.sample = "Stroma")
})
plotting.res <- lapply(Immune.list, function(seurat.x) {
  Subset_marker_plotting(seurat.x, marker.list = Immune.marker.list, CT.sample = "Immune")
})
plotting.res <- lapply(Endometrium.list, function(seurat.x) {
  Subset_marker_plotting(seurat.x, marker.list = Immune.marker.list, CT.sample = "Endometrium")
})
