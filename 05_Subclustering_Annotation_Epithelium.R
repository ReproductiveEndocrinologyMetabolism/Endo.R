#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(pheatmap)

# Setting directories and flags
Project_name = "Endo_All_Epithelium"
Input.dir = "Output/4_Subsetting/"
FindMarker.flag = TRUE
ScaleData.flag = FALSE
SingleR.flag = FALSE
annotation.done = TRUE
theme_set(theme_cowplot())

# Setting the output directory
if (dir.exists(path = paste0("Output/5_Epithelium_Clustering")) == FALSE) {
  print(paste0("Generating output directory Output/5_Epithelium_Clustering"))
  dir.create(path = paste0("Output/5_Epithelium_Clustering"), recursive = TRUE)
  Output.dir = paste0("Output/5_Epithelium_Clustering/")
} else if (dir.exists(path = paste0("Output/5_Epithelium_Clustering")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/5_Epithelium_Clustering/")
} else {
  print("Error with output directory")
}

# Loading reclustered Seurat
print("Seurat object loading")
endo.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_reclustered.h5seurat"))
#endo.integrated = LoadH5Seurat(file = paste0(Output.dir, Project_name, "_reclustered_labelled.h5seurat"))

# Loading immune markers from FindAllMarkers
cluster.markers = readRDS("Output/4_Subsetting/Endo_All_Epithelium_FindAllMarkers.rds")
cluster.markers.top10 = readRDS("Output/4_Subsetting/Endo_All_Epithelium_FindAllMarkers_top10.rds")

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

# Setting markers
### Use the top 3 markers
epithelium_3 = c("EPCAM", "WFDC2", "KLF5", "TACSTD2")
ciliated_3 = c("FOXJ1", "C9orf24", "CDHR3")
Lumenal_epithelium_3.1 = c("PGR", "ESR1", "CPM", "HMGB2", "LGR5")
Lumenal_epithelium_3.2 = c("LGR5", "MMP7", "MKI67", "WNT7A")
Lumenal_epithelium_3.3 = c("LGR5", "SVIL", "DST", "VTCN1")
Glandular_epithelium_3.1 = c("SCGB2A2", "CPM", "ESR1", "HMGB2", "AR")
Glandular_epithelium_3.2 = c("KIAA1324", "SMAD9", "SLC47A1")
SOX9pos = c("MMP7", "FHL2", "MKI67")
SOX9pos_LGRpos = c("WNT7A", "KRT17", "KRT5", "LGR5", "IL32", "PLAU")
SOX9pos_LGRneg = c("IHH", "EMID1", "PPARG", "MMP7", "CPM")
Pre_ciliated = c("MUC12", "CDC20B", "CCNO", "HES6")
Glandular_prolif = c("SLC47A1", "KIAA1324", "SMAD9", "DST", "MMP26")
Lumenal_prolif = c("SVIL", "VTCN1" ,"DST", "IL6")
Early_prolif = c("MMP7", "THBS1", "CADM1")
Late_prolif = c("SCGB1D2", "TRAK1", "ANK3", "ATP1A1")
All_Alonso_Prolif = c("SOX9", "PGR", "ESR1", "MMP7", "CPM", "FHL2", "MKI67", "HMGB2",
                      "PLAU", "IL32", "TNF", "WNT7A", "KRT17", "KRT5", "LGR5",
                      "IHH", "EMID1", "PPARG",
                      "MUC12", "CDC20B", "CCNO", "HES6")
All_Alonso_Secr = c("FOXJ1", "PIFO", "TP73",
                    "HEY1", "ABCH1", "SCGB2A2",
                    "C2CD4A", "SLC18A2", "PAEP", "CXCL14", "SPP1", "GPX3", "DPP4",
                    "PTGS1", "CLDN22", "PAX2", "IL6", "VTCN1", "SLC26A7", "MSLN")
All_Alonso = c("SOX9", "PGR", "ESR1", "MMP7", "CPM", "FHL2", "MKI67", "HMGB2",
               "PLAU", "IL32", "TNF", "WNT7A", "KRT17", "KRT5", "LGR5",
               "IHH", "EMID1", "PPARG",
               "MUC12", "CDC20B", "CCNO", "HES6",
               "FOXJ1", "PIFO", "TP73",
               "HEY1", "ABCH1", "SCGB2A2",
               "C2CD4A", "SLC18A2", "PAEP", "CXCL14", "SPP1", "GPX3", "DPP4",
               "PTGS1", "CLDN22", "PAX2", "IL6", "VTCN1", "SLC26A7", "MSLN")
Figure_curated = c("PTGS1", "CLDN", "PAX", "VTCN1", "SLC26A7",
                   "LGR5", "WNT7A", "KRT5", "SOX9",
                   "CPM", "IHH", "EMID1", "PPARG",
                   "HEY1", "SCGB2A2", "C2CD4A", "SLC18A2", "PAEP", "CXCL14", "SPP1", "GPX3", "DPP4",
                   "MKI67", "HMGB2","AR","CDC20B", "CCNO", 
                   "HES6", "FOXJ1", "PIFO", "TP73")

Epithelium_many = c("SCGB2A2", "CPM", "ESR1", "HMGB2", "AR", "PGR", "EPCAM", "KRT8", "LGR5")
Lumenal_epithelium_many = c("PGR", "ESR1", "CPM", "HMGB2", "LGR5")
Glandular_epithelium_many = c("SCGB2A2", "CPM", "ESR1", "HMGB2", "AR")
Ciliated_epithelium_many = c("MMP7", "HMGB2", "MUC12", "CDC20B", "CCNO", "HES6")
Glandular_secretory_epithelium_many = c("PAEP", "C2CD4A", "SLC18A2", "CXCL14", "GPX3")

# Running markerplotting on seurat clusters with RNA assay
Idents(object = endo.integrated) <- "seurat_clusters"
DefaultAssay(endo.integrated) = "RNA"
Norm.assay = "RNA"

run = Marker_plotting(marker = All_Alonso_Prolif, name = "All_Alonso_Prolif", group = "seurat_clusters")
run = Marker_plotting(marker = All_Alonso_Secr, name = "All_Alonso_Secr", group = "seurat_clusters")

run = Marker_plotting(marker = epithelium_3, name = "epithelium_3", group = "seurat_clusters")
run = Marker_plotting(marker = ciliated_3, name = "ciliated_3", group = "seurat_clusters")
run = Marker_plotting(marker = Lumenal_epithelium_3.1, name = "Lumenal_epithelium_3.1", group = "seurat_clusters")
run = Marker_plotting(marker = Lumenal_epithelium_3.2, name = "Lumenal_epithelium_3.2", group = "seurat_clusters")
run = Marker_plotting(marker = Lumenal_epithelium_3.3, name = "Lumenal_epithelium_3.3", group = "seurat_clusters")
run = Marker_plotting(marker = Glandular_epithelium_3.1, name = "Glandular_epithelium_3.1", group = "seurat_clusters")
run = Marker_plotting(marker = Glandular_epithelium_3.2, name = "Glandular_epithelium_3.2", group = "seurat_clusters")
run = Marker_plotting(marker = SOX9pos, name = "SOX9pos", group = "seurat_clusters")
run = Marker_plotting(marker = SOX9pos_LGRpos, name = "SOX9pos_LGRpos", group = "seurat_clusters")
run = Marker_plotting(marker = SOX9pos_LGRneg, name = "SOX9pos_LGRneg", group = "seurat_clusters")
run = Marker_plotting(marker = Pre_ciliated, name = "Pre_ciliated", group = "seurat_clusters")
run = Marker_plotting(marker = Glandular_prolif, name = "Glandular_prolif", group = "seurat_clusters")
run = Marker_plotting(marker = Lumenal_prolif, name = "Lumenal_prolif", group = "seurat_clusters")
run = Marker_plotting(marker = Early_prolif, name = "Early_prolif", group = "seurat_clusters")
run = Marker_plotting(marker = Late_prolif, name = "Late_prolif", group = "seurat_clusters")

run = Marker_plotting(marker = Epithelium_many, name = "Epithelium_many", group = "seurat_clusters")
run = Marker_plotting(marker = Lumenal_epithelium_many, name = "Lumenal_epithelium_many", group = "seurat_clusters")
run = Marker_plotting(marker = Glandular_epithelium_many, name = "Glandular_epithelium_many", group = "seurat_clusters")
run = Marker_plotting(marker = Ciliated_epithelium_many, name = "Ciliated_epithelium_many", group = "seurat_clusters")
run = Marker_plotting(marker = Glandular_secretory_epithelium_many, name = "Glandular_secretory_epithelium_many", group = "seurat_clusters")

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
  
}

if (annotation.done == TRUE) {
  
  # Labelling with old Epithelium labels
  print("Plotting Epithelium labels")
  endo.integrated$Epithelium_old_labelled = endo.integrated$seurat_clusters
  endo.integrated = SetIdent(endo.integrated, value = "Epithelium_old_labelled")
  endo.integrated = RenameIdents(endo.integrated, "0" = "SOX9 LGR-",
                                 "1" = "Lumenal",
                                 "2" = "Glandular",
                                 "3" = "S-phase proliferative",
                                 "4" = "Lumenal",
                                 "5" = "SOX9 LGR+",
                                 "6" = "SOX9 LGR-",
                                 "7" = "SOX9 LGR-",
                                 "8" = "SOX9 LGR-",
                                 "9" = "SOX9 LGR-",
                                 "10" = "Lumenal",
                                 "11" = "G2M-phase proliferative",
                                 "12" = "SOX9 LGR+",
                                 "13" = "SOX9 LGR-",
                                 "14" = "AR+",
                                 "15" = "Glandular",
                                 "16" = "SOX9 LGR+",
                                 "17" = "SOX9 LGR-",
                                 "18" = "SOX9 LGR+",
                                 "19" = "Ciliated",
                                 "20" = "SOX9 LGR-")
  endo.integrated[["Epithelium_old_labelled"]] = Idents(object = endo.integrated)
  DimPlot(endo.integrated, reduction = "umap", group.by = "Epithelium_old_labelled", 
          label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir, Project_name, "_UMAP_Epithelium_old_labelled.pdf"), dpi = 700)
  
  # Labelling with new Epithelium labels
  print("Plotting Epithelium 230428 labels")
  endo.integrated$Epithelium_labelled_230428 = endo.integrated$seurat_clusters
  endo.integrated = SetIdent(endo.integrated, value = "Epithelium_labelled_230428")
  endo.integrated = RenameIdents(endo.integrated, "0" = "SOX9+ LGR5-",
                                 "1" = "Lumenal",
                                 "2" = "Glandular",
                                 "3" = "SOX9+ proliferative",
                                 "4" = "Lumenal",
                                 "5" = "SOX9+ LGR5+",
                                 "6" = "SOX9+ LGR5-",
                                 "7" = "SOX9+ LGR5-",
                                 "8" = "SOX9+ LGR5-",
                                 "9" = "SOX9+ LGR5-",
                                 "10" = "Lumenal",
                                 "11" = "SOX9+ proliferative",
                                 "12" = "SOX9+ LGR5+",
                                 "13" = "SOX9+ LGR5-",
                                 "14" = "AR+",
                                 "15" = "Glandular",
                                 "16" = "Glandular",
                                 "17" = "SOX9+ LGR5-",
                                 "18" = "SOX9+ LGR5+",
                                 "19" = "Ciliated",
                                 "20" = "SOX9+ LGR5-")
  endo.integrated[["Epithelium_labelled_230428"]] = Idents(object = endo.integrated)
  DimPlot(endo.integrated, reduction = "umap", group.by = "Epithelium_labelled_230428", 
          label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir, Project_name, "_UMAP_Epithelium_labelled_230428.pdf"), dpi = 700)
  
  # Labelling with new Epithelium labels
  print("Plotting Epithelium labels")
  endo.integrated$Epithelium_labelled = endo.integrated$seurat_clusters
  endo.integrated = SetIdent(endo.integrated, value = "Epithelium_labelled")
  endo.integrated = RenameIdents(endo.integrated, "0" = "SOX9+ LGR5-",
                                 "1" = "Lumenal",
                                 "2" = "Glandular",
                                 "3" = "SOX9+ proliferative",
                                 "4" = "SOX9+ LGR5+",
                                 "5" = "SOX9+ LGR5+",
                                 "6" = "SOX9+ LGR5-",
                                 "7" = "SOX9+ LGR5-",
                                 "8" = "SOX9+ LGR5-",
                                 "9" = "SOX9+ LGR5-",
                                 "10" = "SOX9+ LGR5+",
                                 "11" = "SOX9+ proliferative",
                                 "12" = "SOX9+ LGR5+",
                                 "13" = "SOX9+ LGR5-",
                                 "14" = "AR+",
                                 "15" = "Glandular",
                                 "16" = "Lumenal",
                                 "17" = "SOX9+ LGR5-",
                                 "18" = "SOX9+ proliferative",
                                 "19" = "Ciliated",
                                 "20" = "SOX9+ LGR5-")
  endo.integrated[["Epithelium_labelled"]] = Idents(object = endo.integrated)
  DimPlot(endo.integrated, reduction = "umap", group.by = "Epithelium_labelled", 
          label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir, Project_name, "_UMAP_Epithelium_labelled.pdf"), dpi = 700)
  
  DimPlot(endo.integrated, reduction = "umap", split.by = "Group_Stage")
  ggsave2(paste0(Output.dir, Project_name, "_UMAP_Epithelium_labelled_Split-Group.pdf"), dpi = 700,
          height = 12, width = 20)
  
  # Reordering UMAP labels
  Idents(endo.integrated) = "Epithelium_labelled"
  new.order.labels = c("Lumenal", "SOX9+ LGR5+", 
                       "SOX9+ LGR5-", "Glandular", 
                       "SOX9+ proliferative",
                       "AR+", "Ciliated")
  Idents(endo.integrated) <- factor(Idents(endo.integrated), levels= new.order.labels)
  endo.integrated$Epithelium_labelled <- factor(endo.integrated$Epithelium_labelled, levels= new.order.labels)
  
  # 10 colors
  Old.Dimplot.colors = c("darkred", "orangered", 
                     "dodgerblue", "navyblue", 
                     "limegreen", "forestgreen",
                     "deeppink", "darkkhaki")
  
  Dimplot.colors = c("#DD8D61", "#D3020D", 
                     "#2570B7", "#99A3FF", 
                     "#FAA0A1",
                     "#C758C5", "#65C2A5")
  
  # Save the labelled and re-ordered object
  SaveH5Seurat(endo.integrated, paste0(Output.dir,Project_name, "_reclustered_labelled.h5seurat"), overwrite = TRUE)
  
  ##Idents(endo.copy) <- factor(Idents(endo.copy), levels= new.order.labels)
  ###endo.copy$Immune_labelled <- factor(endo.copy$Immune_labelled, levels= new.order.labels)
  ###Dimplot.color = RColorBrewer::brewer.pal(16, "Set1")
  DimPlot(endo.integrated, cols = Dimplot.colors, raster = FALSE)
  ggsave2(paste0(Output.dir, Project_name, "_no-labels_UMAP_undefined_figure_ptdefault.pdf"), dpi = 700)
  
  DimPlot(endo.integrated, cols = Dimplot.colors, label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir, Project_name, "_labels_UMAP_undefined_figure.pdf"), dpi = 700)
  
  
  DotPlot(object = endo.integrated, features = All_Alonso, group.by = "Epithelium_labelled",
          cols = c("dodgerblue", "firebrick"), dot.min = 0.1, assay = "RNA", dot.scale = 10) + 
    theme(axis.text.x=element_text(angle=90, vjust = 0.3, hjust = 1))
  ggsave2(paste0(Output.dir, Project_name, "_Dotplot_Epithelium_Curated.pdf"), dpi = 700, width = 12)
  
  SaveH5Seurat(endo.integrated, paste0(Output.dir, Project_name, "_reclustered_labelled.h5seurat"), overwrite = TRUE)
  
  DotPlot(object = endo.integrated, features = All_Alonso, group.by = "Epithelium_labelled",
          cols = c("navajowhite", "firebrick"), dot.min = 0, dot.scale = 10) + 
    theme(axis.text.x=element_text(angle=90))
  ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Marker_Dotplot.pdf"), 
          dpi = 700, width = 14)
  
  endo.average = AverageExpression(endo.integrated, return.seurat = TRUE)
  
  DoHeatmap(endo.average, features = All_Alonso, raster = FALSE,
            group.bar = TRUE, group.colors = c(rep("white", length(levels(endo.average)))),
            draw.lines = FALSE, angle = 45, hjust = 0) + 
    scale_fill_gradientn(colors = c("dodgerblue", "seashell", "tomato", "firebrick"))
  
  ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Marker_Heatmap.pdf"), 
          dpi = 700, height = 14, width = 8)
  
  DoHeatmap(endo.average, features = Figure_curated, raster = FALSE,
            group.bar = TRUE, group.colors = c(rep("white", length(levels(endo.average)))),
            draw.lines = FALSE, angle = 45, hjust = 0) + 
    scale_fill_gradientn(colors = c("dodgerblue", "seashell", "tomato", "firebrick"))
  
  ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Marker_Heatmap_selected_genes.pdf"), 
          dpi = 700, height = 14, width = 8)
  
  # Plotting the integrated object
  print("Plotting of UMAP")
  UMAP_Ident = DimPlot(endo.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
  ggsave2(paste0(plot = UMAP_Ident, Output.dir,Project_name,"_labelled_UMAP_Ident.pdf"), dpi = 700)
  UMAP_Phase = DimPlot(endo.integrated, reduction = "umap", group.by = "Phase")
  ggsave2(paste0(plot = UMAP_Phase, Output.dir,Project_name,"_UMAP_Phase.pdf"), dpi = 700)
  UMAP_Label = DimPlot(endo.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
  ggsave2(plot = UMAP_Label, paste0(Output.dir,Project_name,"_labelled_UMAP_Clusters.pdf"), dpi = 700)
  UMAP_Group = DimPlot(endo.integrated, reduction = "umap", group.by = "Group_Stage",)
  ggsave2(plot = UMAP_Group, paste0(Output.dir,Project_name,"_UMAP_Group_Stage.pdf"), dpi = 700)
  print("Done with basic UMAP")
}


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
Plot_CQ_Epithelium = Plotting_QC(x = endo.integrated, group = "Epithelium_labelled")

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

if (FindMarker.flag == TRUE) {
  # Setting cluster to be tested to SCT clustered
  Idents(object = endo.integrated) <- "Epithelium_labelled"
  DefaultAssay(endo.integrated) = "RNA"
  Norm.assay = "RNA"
  
  print("Finding markers")
  endo.markers <-FindAllMarkers(endo.integrated, assay = "RNA")
  print("Saving markers as .csv")
  write.csv(endo.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_FindAllMarkers.csv"), quote = F)
  print("Saving markers as .rds")
  saveRDS(endo.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_Endo_markers.rds"))
  print("Extracting top 10 markers per cluster")
  top10 = endo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  
  print("Generating heatmap")
  endo.heatmap = DoHeatmap(endo.integrated, features = top10$gene) + NoLegend()
  ggsave2(filename = paste0(Output.dir, Project_name, "_", Norm.assay, "_Top10_Genes_Heatmap.pdf"),
          plot = endo.heatmap,
          dpi = 700)
  
  
}


