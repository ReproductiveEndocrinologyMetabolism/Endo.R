#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(pheatmap)

# Setting directories and flags
Project_name = "Endo_All_Immune"
Input.dir = "Output/4_Subsetting/"
FindMarker.flag = FALSE
ScaleData.flag = FALSE
SingleR.flag = TRUE
theme_set(theme_cowplot())

# Setting the output directory
if (dir.exists(path = paste0("Output/5_Immune_Clustering")) == FALSE) {
  print(paste0("Generating output directory Output/5_Immune_Clustering"))
  dir.create(path = paste0("Output/5_Immune_Clustering"), recursive = TRUE)
  Output.dir = paste0("Output/5_Immune_Clustering/")
} else if (dir.exists(path = paste0("Output/5_Immune_Clustering")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/5_Immune_Clustering/")
} else {
  print("Error with output directory")
}

# Loading reclustered Seurat
print("Seurat object loading")
#endo.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_reclustered.h5seurat"))
endo.integrated = LoadH5Seurat(file = paste0(Output.dir, Project_name, "_reclustered_labelled.h5seurat"))

# Loading immune markers from FindAllMarkers
immune.markers = readRDS("Output/4_Subsetting/Endo_All_Immune_FindAllMarkers.rds")
immune.markers.top10 = readRDS("Output/4_Subsetting/Endo_All_Immune_FindAllMarkers_top10.rds")

# Data is re-scaled after subsetting as the mean and SD will have changed 
if (ScaleData.flag == TRUE) {
  print("Setting DefaultAssay to RNA and log2 normalising it")
  DefaultAssay(endo.integrated) = "RNA"
  endo.integrated = NormalizeData(endo.integrated)
  endo.integrated = FindVariableFeatures(endo.integrated)
  all.genes = rownames(endo.integrated)
  endo.integrated = ScaleData(endo.integrated, features = all.genes)  #endo.integrated = FindVariableFeatures(endo.integrated)
  
} else if (ScaleData.flag == FALSE) {
  print("Setting DefaultAssay to RNA")
  DefaultAssay(endo.integrated) = "RNA" #Alternative is RNA
}

# Plotting function
Marker_plotting <- function(marker, name, group, selected.assay = Norm.assay) {
  
  print(paste("Plotting", name, "with", selected.assay))
  
  DefaultAssay(endo.integrated) = selected.assay
  VlnPlot(endo.integrated, features = marker, pt.size = 0, sort = "increasing", group.by = group)
  ggsave2(paste0(Output.dir, Project_name, "_", group, "_", name, "_", Norm.assay, "_Marker_ViolinPlot.pdf"), dpi = 700)
  
  FeaturePlot(endo.integrated, features = marker)
  ggsave2(paste0(Output.dir, Project_name, "_", group, "_", name, "_", Norm.assay, "_Marker_FeaturePlot.pdf"), dpi = 700)
  
  DotPlot(object = endo.integrated, features = marker, group.by = group) + RotatedAxis()
  ggsave2(paste0(Output.dir, Project_name, "_", group, "_", name, "_", Norm.assay, "_Marker_Dotplot.pdf"), dpi = 700)
}

# Curated gene marker list:
uNK_1_curated = c("NCAM1", "SPINK2", "CSF1", "CYP26A1", "B4GALNT1", "CSF1", 
                  "GZMA", "EPAS1", "STAT3", "GNLY", "PRF1")
uNK_2_curated = c("NCAM1", "CD2", "ITGB2", "ANXA1", "CD27", "GZMH")
uNK_3_curated = c("NCAM1", "CD160", "CCL5", "CXCR4", "TIGIT")
Peripheral_NK_curated = c("FCGR3A", "CD160")
uNK_1_curated = c("NCAM1", "HMGB2", "MKI67")
uNK_2_curated = c("NCAM1", "CD2", "ITGB2", "ANXA1", "CD27", "GZMH")
uNK_3_curated = c("NCAM1", "CD160", "CCL5")
uNK_all_curated = c("NCAM1", "HMGB2", "MKI67", "ITGA1", "SPINK2", "CSF1", "CYP26A1", "B4GALNT1",
                    "CD2", "ITGB2", "ANXA1", "CD27", "GZMH",
                    "CD160", "CCL5", "TIGIT", "CD69","CXCR4",
                    "FCGR3A")
cycling.lymphocytes_ILC3_pNK_uNK_mast.cells = c("CD69", "CD3G", "CD8A", "MKI67", "HMGB2",
                                                "MS4A2", "KIT", "IL7R", "NCAM1", "FCGR3A",
                                                "CXCR4", "CD160", "CCL5", "KIR2DL1", "SPINK2",
                                                "CSF1")
uM_curated = c("CD14", "CD68", "EREG", "IL1B", "SELENOP", "HMOX1")
Tcells_CD4_curated = c("CD69", "CD3G","IL7R", "CCR5", "CD4",
                       "GZMK", "IL10", "GZMA", "SOST", "IL22RA2",
                       "SIGLEC17P", "IL26", "PENK", "MIR144")
Tcells_CD8_curated = c("CD69", "CD3G","CCL5", "PSMB10", "GZMA", "GZMH", "PRF1", "CD8A") #Yan, Lu, Yan, Wang Frontiers in Immunology 2021
Tcells_CD4_CD8_curated = c("CD69", "CD3G", "IL7R", "CCR5", "CD4",
                           "CCL5", "PSMB10", "GZMA", "GZMH", "PRF1", "CD8A",
                           "CD27", "CCR7", "IKZF4", #Treg
                           "ZNF683", "GNG4", "PDCD1", #CD8 a/a
                           "TOX2", "SATB1", "CCR9", # CD8 a/b
                           "GZMK", "IL10") # CD4 

DC_1_curated = c("CD14", "CD68", "IL1B", "BATF3", "CADM1", "CLEC9A")
DC_2_curated = c("CD14", "CD68", "SELENOP", "CLEC10A", "FCER1A", "CD1C")
Cycling_DC_curated = c("MKI67", "TOP2A", "CLEC10A")
Migratory_DC_curated = c("EBI3", "CCR7", "CCL19")
pDC_curated = c("IL3RA", "LILRA4", "PLD4")
DCs = c("CD14", "CSF1R", "CD1C", "CD68", "IL1B", "HMGB2", "CXCR4")
uterine_Macrophages = c("CD69", "CD14", "EREG", "CSF1R", "CD68", "HMOX1", "IL1B", "SELENOP", "FCGR3A")

DC_curated = c("CD1C", "IL1B", "CD14", "CD68", "BATF3", "CADM1", "CLEC9A",
               "CLEC10A", "FCER1A",
               "MKI67", "TOP2A",
               "EBI3", "CCR7", "CCL19",
               "IL3RA", "LILRA4", "PLD4",
               "CD1A", "CD83",
               "IRF8", "PCLAF") #precursor DC

B.cells_curated = c("CD79A", "MS4A1", "CD19",
                    "CXCR5", "TNFRSF13B", "CD22",
                    "IGHM", "IGHD")
Mast_cells_curated = c("HPGDS", "TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT") 
ILC3_curated = c("IL7R", "S100A13", "TLE1", "AREG",
                 "CXCR3", "CD3D", "IKZF3",
                 "GATA3", "KLRG1", "HPGDS", 
                 "IL4I1", "RORC", "KIT")
ILCs_curated = c("S100A13", "TLE1", "AREG",
                 "CXCR3", "CD3D", "IKZF3",
                 "GATA3", "KLRG1", "HPGDS", 
                 "IL4I1", "RORC", "KIT")

# All curated markers in the order of DC1, DC2, cycling DC, migratory DC,
# pDC, uM1, uM2, T-cell CD8, T-cell CD4, B-cell, uNK, mast cells
Immune_all_curated = c("BATF3", "CADM1", "CLEC9A", #DC1
                       "CLEC10A", "FCER1A", "CD1C", #DC2
                       "MKI67", "TOP2A", #Cycling DC
                       "EBI3", "CCR7", "CCL19", # Migratory DC
                       "IL3RA", "LILRA4", "PLD4", #pDC
                       "CD14", "CD68", "IL1B",  #uM1
                       "SELENOP", "HMOX1", #uM2
                       "CD69", "CD3G", "CCL5", "CD8A", "IL7R", "CD4", #T-cells
                       "ZCCHC7", "MS4A1", "IGHM", # B-cells
                       "NCAM1", "ITGA1", "GNLY", # uNK
                       "CPA3", "MS4A2", "KIT") #Mast cells

# Plotting heatmap and dotplot of all curated markers
Immune_all_curated_figure = c("FOXP3", "CD3G", "CD27", "IKZF4", "CD8A", "CCL5", "CD4", 
                              "IL7R", "RORC", #Treg, T-cell CD8, T-cell CD4, ILC3
                              "MS4A1", "IGHM", #B-celler
                              "NCAM1", "ITGA1", "SPINK2", "CSF1", "CD160", "GNLY", # uNK 1-3
                              "IL3RA", "LILRA4", "PLD4",# pDC
                              "EBI3", "CCR7", "CCL19", # Migratory DC
                              "BATF3", "CADM1", "CLEC9A", #DC1
                              "CLEC10A", "FCER1A", "CD1C", #DC2
                              "CD14", "SELENOP", "HMOX1", "IL1B", #uM 1-2
                              "CPA3", "KIT", "MS4A2") # Mast cells

# Loading top markers based on FindAllMarkers per seurat cluster
immune.markers.labelled = readRDS("Output/5_Immune_Clustering/Endo_All_Immune_Endo_markers_Immune_labelled.rds")

# Running markerplotting on seurat clusters with RNA assay
Idents(object = endo.integrated) <- "seurat_clusters"
DefaultAssay(endo.integrated) = "RNA"
Norm.assay = "RNA"

run = Marker_plotting(marker = Immune, name = "Immune_clusters", group = "seurat_clusters")
run = Marker_plotting(marker = Lymphoid, name = "Lymphoid_clusters", group = "seurat_clusters")
run = Marker_plotting(marker = Myehoid, name = "Myehoid_clusters", group = "seurat_clusters")

# Plotting only curated markers
DefaultAssay(endo.integrated) = "RNA"
Norm.assay = "RNA"
run = Marker_plotting(marker = uNK_all_curated, name = "uNK_all_curated", group = "seurat_clusters")
run = Marker_plotting(marker = Vento_Tormo_uNK, name = "Vento_Tormo_uNK", group = "seurat_clusters")
run = Marker_plotting(marker = Vento_Tormo_uNK_supp, name = "Vento_Tormo_uNK_supp", group = "seurat_clusters")
run = Marker_plotting(marker = uM_curated, name = "uM_curated", group = "seurat_clusters")
run = Marker_plotting(marker = Tcells_CD4_CD8_curated, name = "Tcells_CD4_CD8_curated", group = "seurat_clusters")
run = Marker_plotting(marker = DC_curated, name = "DC_curated", group = "seurat_clusters")
run = Marker_plotting(marker = B.cells_curated, name = "B.cells_curated", group = "seurat_clusters")
run = Marker_plotting(marker = Immune_all_curated, name = "Immune_all_curated", group = "seurat_clusters")
run = Marker_plotting(marker = Mast_cells_curated, name = "Mast_cells_curated", group = "seurat_clusters")
run = Marker_plotting(marker = ILC3_curated, name = "ILC3_curated", group = "seurat_clusters")
run = Marker_plotting(marker = ILCs_curated, name = "ILC-all_curated", group = "seurat_clusters")
run = Marker_plotting(marker = cycling.lymphocytes_ILC3_pNK_uNK_mast.cells, 
                      name = "cycling-lymphocytes_ILC3_pNK_uNK_mast-cells", group = "seurat_clusters")

# Do automated annotation
if (SingleR.flag == TRUE) {
  
  library(SingleR)
  library(celldex)
  
  # Setting annotation for SingleR annotation
  SingleR.annotation.ref.Monaco = celldex::MonacoImmuneData()
  SingleR.annotation.ref.DbImmune = celldex::DatabaseImmuneCellExpressionData(cell.ont = "nonna")
  SingleR.annotation.ref.Blueprint = celldex::BlueprintEncodeData(cell.ont = "nonna")
  
  endo.sce = as.SingleCellExperiment(DietSeurat(endo.integrated))
  
  # Run automated annotation using several different references
  singleR.main.Monaco = SingleR(test = endo.sce, assay.type.test = 1, ref = SingleR.annotation.ref.Monaco,labels = SingleR.annotation.ref.Monaco$label.main)
  singleR.main.DbImmune = SingleR(test = endo.sce, assay.type.test = 1, ref = SingleR.annotation.ref.DbImmune,labels = SingleR.annotation.ref.DbImmune$label.main)
  singleR.main.Blueprint = SingleR(test = endo.sce, assay.type.test = 1, ref = SingleR.annotation.ref.Blueprint,labels = SingleR.annotation.ref.Blueprint$label.main)
  
  singleR.fine.Monaco = SingleR(test = endo.sce, assay.type.test = 1, ref = SingleR.annotation.ref.Monaco,labels = SingleR.annotation.ref.Monaco$label.fine)
  singleR.fine.DbImmune = SingleR(test = endo.sce, assay.type.test = 1, ref = SingleR.annotation.ref.DbImmune,labels = SingleR.annotation.ref.DbImmune$label.fine)
  singleR.fine.Blueprint = SingleR(test = endo.sce, assay.type.test = 1, ref = SingleR.annotation.ref.Blueprint,labels = SingleR.annotation.ref.Blueprint$label.fine)
  

  endo.SingleR = endo.integrated
  
  # Testing and plotting Monaco main
  endo.SingleR$immune_labels_monaco_main = singleR.main.Monaco$pruned.labels
  Idents(endo.SingleR) = "immune_labels_monaco_main"
  DimPlot(endo.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_Immune_cells_SingleR_Monaco_main.pdf"), dpi = 700)
  
  # Testing and plotting DbImmune main
  endo.SingleR$immune_labels_DbImmune_main = singleR.main.DbImmune$pruned.labels
  Idents(endo.SingleR) = "immune_labels_DbImmune_main"
  DimPlot(endo.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_Immune_cells_SingleR_DbImmune_main.pdf"), dpi = 700)
  
  # Testing and plotting Blueprint main
  endo.SingleR$immune_labels_Blueprint_main = singleR.main.Blueprint$pruned.labels
  Idents(endo.SingleR) = "immune_labels_Blueprint_main"
  DimPlot(endo.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_Immune_cells_SingleR_Blueprint_main.pdf"), dpi = 700)
  
  # Testing and plotting Monaco fine
  endo.SingleR$immune_labels_monaco_fine = singleR.fine.Monaco$pruned.labels
  Idents(endo.SingleR) = "immune_labels_monaco_fine"
  DimPlot(endo.SingleR, label = T)
  DimPlot(endo.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_Immune_cells_SingleR_Monaco_fine.pdf"), dpi = 700,
          width = 20, height = 20)
  
  # Testing and plotting DbImmune fine
  endo.SingleR$immune_labels_DbImmune_fine = singleR.fine.DbImmune$pruned.labels
  Idents(endo.SingleR) = "immune_labels_DbImmune_fine"
  DimPlot(endo.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_Immune_cells_SingleR_DbImmune_fine.pdf"), dpi = 700)
  
  # Testing and plotting Blueprint fine
  endo.SingleR$immune_labels_Blueprint_fine = singleR.fine.Blueprint$pruned.labels
  Idents(endo.SingleR) = "immune_labels_Blueprint_fine"
  DimPlot(endo.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_Immune_cells_SingleR_Blueprint_fine.pdf"), dpi = 700)
  
  
  SingleR.annotation.ref.Monaco = celldex::MonacoImmuneData()
  SingleR.annotation.ref.DbImmune = celldex::DatabaseImmuneCellExpressionData("nonna")
  SingleR.annotation.ref.Blueprint = celldex::BlueprintEncodeData("nonna")
  
  # Running markerplotting on SingleR annotation and RNA
  #Idents(object = endo.integrated) <- "monaco.main"
  DefaultAssay(endo.integrated) = "RNA"
  Norm.assay = "RNA"
  
  ######### PLOTTING SINGLE R MARKERS #########
  run = Marker_plotting(marker = Immune, name = "Immune_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = Lymphoid, name = "Lymphoid_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = Myehoid, name = "Myehoid_SingleR", group = "immune_labels")
  
  # Lymphoids:
  run = Marker_plotting(marker = uNK, name = "Lymphoid_uNK_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = Peripheral_NK, name = "Lymphoid_Peripheral_NK_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = B_cells, name = "Lymphoid_B_cells_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = B_cells_3, name = "Lymphoid_B_cells_3_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = Tcells_CD8, name = "Lymphoid_Tcells_CD8_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = Tcells_CD4, name = "Lymphoid_Tcells_CD4_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = NKT, name = "NKT", group = "immune_labels")
  run = Marker_plotting(marker = Innate_lymphoid_cells, name = "Lymphoid_Innate_lymphoid_cells_SingleR", group = "immune_labels")
  # Myeloids: 
  run = Marker_plotting(marker = uterine_Macrophages, name = "Myehoid_uterine_Macrophages_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = mast_cells, name = "Myehoid_mast_cells_SingleR", group = "immune_labels")
  # Other:
  run = Marker_plotting(marker = DCs, name = "DCs_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = cycling_lymphocytes, name = "cycling_lymphocytes_SingleR", group = "immune_labels")
  run = Marker_plotting(marker = All_immune, name = "All_immune_SingleR", group = "immune_labels")
  
  #endo.integrated <- SetIdent(endo.integrated, value = "SingleR.annotation.fine")
  #DimPlot(endo.integrated, label = T , repel = T, label.size = 3) + NoLegend()
  ########### END OF PLOTTING SINGLE R MARKERS ############
}

# Labelling with new immune labels
print("Plotting Immune labels")
endo.integrated$Immune_labelled = endo.integrated$seurat_clusters
endo.integrated = SetIdent(endo.integrated, value = "Immune_labelled")
endo.integrated = RenameIdents(endo.integrated, "0" = "uNK 2",
                               "1" = "uNK 1",
                               "2" = "T-cells CD8+",
                               "3" = "uM 1",
                               "4" = "T-cells CD4+",
                               "5" = "uNK 3",
                               "6" = "uM 2",
                               "7" = "uM 2",
                               "8" = "Undefined",
                               "9" = "uM 1",
                               "10" = "uNK 2",
                               "11" = "uM 1",
                               "12" = "DC2",
                               "13" = "Tregs",
                               "14" = "DC1",
                               "15" = "Undefined",
                               "16" = "uNK 1",
                               "17" = "Mast cells",
                               "18" = "B-cells",
                               "19" = "uNK 1",
                               "20" = "ILC3",
                               "21" = "Migratory DC",
                               "22" = "pDC")
endo.integrated[["Immune_labelled"]] = Idents(object = endo.integrated)
DimPlot(endo.integrated, reduction = "umap", group.by = "Immune_labelled", 
        label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir, Project_name, "_UMAP_Immune_labelled.pdf"), dpi = 700)

DimPlot(endo.integrated, reduction = "umap", split.by = "Group_Stage")
ggsave2(paste0(Output.dir, Project_name, "_UMAP_Immune_labelled_Split-Group.pdf"), dpi = 700,
        height = 12, width = 20)

# Reordering UMAP labels
Idents(endo.integrated) = "Immune_labelled"
new.order.labels = c("T-cells CD4+", "T-cells CD8+", 
                     "uNK 1", "uNK 2", "uNK 3",
                     "uM 1", "uM 2", "DC1", "DC2",
                     "Tregs", "ILC3", "B-cells", 
                     "Mast cells", 
                     "Migratory DC", "pDC", "Undefined")
Idents(endo.integrated) <- factor(Idents(endo.integrated), levels= new.order.labels)
endo.integrated$Immune_labelled <- factor(endo.integrated$Immune_labelled, levels= new.order.labels)

#Dimplot.colors = c("darkred", "orangered", # T-cells
#                   "dodgerblue", "skyblue", "navyblue", # uNKs
#                   "sandybrown", "darkkhaki", # uMs
#                   "salmon", #Treg
#                   "violetred", # ILC3
#                   "deeppink", "maroon", # B-cells, mast cells
#                   "seagreen", "limegreen", #mDC
#                   "thistle", "rosybrown", # migratory DC, pDC
#                   "dimgrey") # Undefined

Dimplot.colors = c("darkred", "orangered", # T-cells
                   "dodgerblue", "skyblue", "navyblue", # uNKs
                   "plum", "maroon", # uMs
                   "seagreen", "limegreen", #mDC
                   "salmon", #Treg
                   "violetred", # ILC3
                   "sandybrown", "darkkhaki", # B-cells, mast cells
                   "palegoldenrod", "rosybrown", # migratory DC, pDC
                   "dimgrey") # Undefined

# Save the labelled and re-ordered object
SaveH5Seurat(endo.integrated, paste0(Output.dir,Project_name, "_reclustered_labelled.h5seurat"), overwrite = TRUE)

##Idents(endo.copy) <- factor(Idents(endo.copy), levels= new.order.labels)
###endo.copy$Immune_labelled <- factor(endo.copy$Immune_labelled, levels= new.order.labels)
###Dimplot.color = RColorBrewer::brewer.pal(16, "Set1")
DimPlot(endo.integrated, cols = Dimplot.colors, raster = FALSE)
ggsave2(paste0(Output.dir, Project_name, "_no-labels_UMAP_undefined_figure_ptdefault.pdf"), dpi = 700)

DimPlot(endo.integrated, cols = Dimplot.colors, raster = FALSE, pt.size = 0.5)
ggsave2(paste0(Output.dir, Project_name, "_no-labels_UMAP_undefined_figure_pt0.5.pdf"), dpi = 700)

DimPlot(endo.integrated, cols = Dimplot.colors, label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir, Project_name, "_labels_UMAP_undefined_figure.pdf"), dpi = 700)


DotPlot(object = endo.integrated, features = Immune_all_curated, group.by = "Immune_labelled",
        cols = c("blue", "red"), dot.min = 0.05, assay = "RNA", dot.scale = 10) + 
  theme(axis.text.x=element_text(angle=90, vjust = 0.3, hjust = 1))
ggsave2(paste0(Output.dir, Project_name, "_Dotplot_Immune_Curated.pdf"), dpi = 700, width = 12)

SaveH5Seurat(endo.integrated, paste0(Output.dir, Project_name, "_reclustered_labelled.h5seurat"), overwrite = TRUE)

# Running markerplotting on Immune clusters with SRNA assay
Idents(object = endo.integrated) <- "Immune_labelled"
DefaultAssay(endo.integrated) = "RNA"
Norm.assay = "RNA"

# Plotting only curated markers
run = Marker_plotting(marker = uNK_all_curated, name = "uNK_all_curated", group = "Immune_labelled")
run = Marker_plotting(marker = Vento_Tormo_uNK, name = "Vento_Tormo_uNK", group = "Immune_labelled")
run = Marker_plotting(marker = Vento_Tormo_uNK_supp, name = "Vento_Tormo_uNK_supp", group = "Immune_labelled")
run = Marker_plotting(marker = uM_curated, name = "uM_curated", group = "Immune_labelled")
run = Marker_plotting(marker = Tcells_CD4_CD8_curated, name = "Tcells_CD4_CD8_curated", group = "Immune_labelled")
run = Marker_plotting(marker = DC_curated, name = "DC_curated", group = "Immune_labelled")
run = Marker_plotting(marker = B.cells_curated, name = "B.cells_curated", group = "Immune_labelled")
run = Marker_plotting(marker = Immune_all_curated, name = "Immune_all_curated", group = "Immune_labelled")
run = Marker_plotting(marker = Mast_cells_curated, name = "Mast_cells_curated", group = "Immune_labelled")
run = Marker_plotting(marker = ILC3_curated, name = "ILC3_curated", group = "Immune_labelled")
run = Marker_plotting(marker = ILCs_curated, name = "ILC-all_curated", group = "Immune_labelled")
run = Marker_plotting(marker = cycling.lymphocytes_ILC3_pNK_uNK_mast.cells, 
                      name = "cycling-lymphocytes_ILC3_pNK_uNK_mast-cells", group = "Immune_labelled")

DotPlot(object = endo.integrated, features = Immune_all_curated, group.by = "Immune_labelled",
        cols = c("blue", "red"), dot.min = 0.1) + theme(axis.text.x=element_text(angle=90))
ggsave2(paste0(Output.dir, Project_name, "_", group, "_", name, "_", Norm.assay, "_Marker_Dotplot.pdf"), dpi = 700)

run = Marker_plotting(marker = Immune_all_curated_figure, name = "Immune_all_curated", group = "Immune_labelled")

gene.marker.plot.order = rev(c("Tregs", "T-cells CD8+", 
                               "T-cells CD4+", "ILC3", "B-cells",
                               "uNK 1", "uNK 2", "uNK 3",
                               "pDC", "Migratory DC", 
                               "DC1", "DC2", "uM 1", "uM 2",
                               "Mast cells", "Undefined"))
Idents(endo.integrated) <- factor(Idents(endo.integrated), levels= gene.marker.plot.order)
endo.integrated$Immune_labelled <- factor(endo.integrated$Immune_labelled, levels= gene.marker.plot.order)

DotPlot(object = endo.integrated, features = Immune_all_curated_figure, group.by = "Immune_labelled",
       cols = c("navajowhite", "firebrick"), dot.min = 0, dot.scale = 10) + 
  theme(axis.text.x=element_text(angle=90))
ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Marker_Dotplot.pdf"), 
        dpi = 700, width = 14)


gene.marker.plot.order = c("Tregs", "T-cells CD8+", 
                           "T-cells CD4+", "ILC3", "B-cells",
                           "uNK 1", "uNK 2", "uNK 3",
                           "pDC", "Migratory DC", 
                           "DC1", "DC2", "uM 1", "uM 2",
                           "Mast cells", "Undefined")


Idents(endo.integrated) <- factor(Idents(endo.integrated), levels= gene.marker.plot.order)
endo.integrated$Immune_labelled <- factor(endo.integrated$Immune_labelled, levels= gene.marker.plot.order)

endo.average = AverageExpression(endo.integrated, return.seurat = TRUE)

DoHeatmap(endo.average, features = Immune_all_curated_figure, draw.lines = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) 

DoHeatmap(endo.average, features = Immune_all_curated_figure, raster = FALSE,
          group.bar = TRUE, group.colors = c(rep("white", length(levels(endo.average)))),
          draw.lines = FALSE, angle = 45, hjust = 0) +
  theme(axis.text.y = element_text(face = "italic")) +    
  scale_fill_gradientn(colors = c("dodgerblue", "seashell", "tomato", "firebrick"))


ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Marker_Heatmap.pdf"), 
        dpi = 700, height = 14, width = 8)

# Plotting the integrated object
print("Plotting of UMAP")
UMAP_Ident = DimPlot(endo.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir,Project_name,"_labelled_UMAP_Ident.pdf"), dpi = 700)
UMAP_Phase = DimPlot(endo.integrated, reduction = "umap", group.by = "Phase")
ggsave2(paste0(Output.dir,Project_name,"_UMAP_Phase.pdf"), dpi = 700)
UMAP_Label = DimPlot(endo.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir,Project_name,"_labelled_UMAP_Clusters.pdf"), dpi = 700)
UMAP_Group = DimPlot(endo.integrated, reduction = "umap", group.by = "Group_Stage",)
ggsave2(paste0(Output.dir,Project_name,"_UMAP_Group_Stage.pdf"), dpi = 700)
UMAP_Group = DimPlot(endo.integrated, reduction = "umap", group.by = "Group_Treatment")
ggsave2(paste0(Output.dir,Project_name,"_UMAP_Group_Treatment.pdf"), dpi = 700)
print("Done with basic UMAP")

# Plotting QC metric between groups, samples and clusters
Plotting_QC <- function(x, group = "") {
  
  # Generating ridgeplots
  RidgePlot(endo.integrated, features = "percent.mt", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_",group, "_Ridgeplot_mtDNA.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "nCount_RNA", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_",group, "_Ridgeplot_nCount_RNA.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "nFeature_RNA", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_",group, "_Ridgeplot_nFeature_RNA.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "S.Score", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_",group, "_Ridgeplot_S-Score.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "G2M.Score", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_",group, "_Ridgeplot_G2M-Score.pdf"), dpi = 700)
  
  # Generating violin plots
  VlnPlot(endo.integrated, features = "percent.mt", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_",group, "_Vlnplot_mtDNA.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "nCount_RNA", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_",group, "_Vlnplot_nCount_RNA.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "nFeature_RNA", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_",group, "_Vlnplot_nFeature_RNA.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "S.Score", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_",group, "_Vlnplot_S-Score.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "G2M.Score", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_",group, "_Vlnplot_G2M-Score.pdf"), dpi = 700)
  
}

Plot_CQ_Group.Stage = Plotting_QC(x = endo.integrated, group = "Group_Stage")
Plot_CQ_Group.Stage = Plotting_QC(x = endo.integrated, group = "Group_Treatment")
Plot_CQ_Pat_nr = Plotting_QC(x = endo.integrated, group = "orig.ident")
Plot_CQ_seurat_clusters = Plotting_QC(x = endo.integrated, group = "seurat_clusters")
Plot_CQ_seurat_clusters = Plotting_QC(x = endo.integrated, group = "Immune_labelled")

# Generating feature plots
FeaturePlot(endo.integrated, features = "percent.mt")
ggsave2(paste0(Output.dir,Project_name, "_FeaturePlot_mtDNA.pdf"), dpi = 700)

FeaturePlot(endo.integrated, features = "nCount_RNA")
ggsave2(paste0(Output.dir,Project_name, "_FeaturePlot_nCount_RNA.pdf"), dpi = 700)

FeaturePlot(endo.integrated, features = "nFeature_RNA")
ggsave2(paste0(Output.dir,Project_name, "_FeaturePlot_nFeature_RNA.pdf"), dpi = 700)

FeaturePlot(endo.integrated, features = "S.Score")
ggsave2(paste0(Output.dir,Project_name, "_FeaturePlot_S-Score.pdf"), dpi = 700)

FeaturePlot(endo.integrated, features = "G2M.Score")
ggsave2(paste0(Output.dir,Project_name, "_FeaturePlot_G2M-Score.pdf"), dpi = 700)