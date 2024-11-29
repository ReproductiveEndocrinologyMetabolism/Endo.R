#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(pheatmap)

# Setting directories and flags
Project_name = "Endo_All_Endothelial"
Input.dir = "Output/4_Subsetting/"
FindMarker.flag = TRUE
ScaleData.flag = FALSE
SingleR.flag = FALSE
annotation.done = TRUE
theme_set(theme_cowplot())

# Setting the output directory
if (dir.exists(path = paste0("Output/5_Endothelial_Clustering")) == FALSE) {
  print(paste0("Generating output directory Output/5_Endothelial_Clustering"))
  dir.create(path = paste0("Output/5_Endothelial_Clustering"), recursive = TRUE)
  Output.dir = paste0("Output/5_Endothelial_Clustering/")
} else if (dir.exists(path = paste0("Output/5_Endothelial_Clustering")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/5_Endothelial_Clustering/")
} else {
  print("Error with output directory")
}

# Loading reclustered Seurat
print("Seurat object loading")
endo.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_reclustered.h5seurat"))
#endo.integrated = LoadH5Seurat(file = paste0(Output.dir, Project_name, "_reclustered_labelled.h5seurat"))

# Loading immune markers from FindAllMarkers
cluster.markers = readRDS("Output/4_Subsetting/Endo_All_Endothelial_FindAllMarkers.rds")
cluster.markers.top10 = readRDS("Output/4_Subsetting/Endo_All_Endothelial_FindAllMarkers_top10.rds")

UMAP_original_label = DimPlot(endo.integrated, reduction = "umap", group.by = "Labelled_Clusters_SCT.1",)
ggsave2(plot =  UMAP_original_label, filename = paste0(Output.dir,Project_name,"_Labelled_Clusters_SCT.1.pdf"), dpi = 700)
print("Done with basic UMAP")  

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

# Setting markers
### Use the top 3 markers
Alonso_Endothelial_All = c("ACTA2", "MYH11", "GUCY1A2", "IGF1", 
                           "CD34", "SEMA3G", "GJA5", 
                           "ACKR1", "PLVAP",
                           "PROX1", "FLT4")
Alonso_Artery_vein = c("CD34", "SEMA3G", "GJA5", "ACKR1", "PLVAP")
VEGF_genes = c("VEGFA", "VEGFB")
VEGF_receptor = c("FLT1", "KDR", "FLT4")
Hormones = c("AR", "ESR1", "PGR")
Horomones_VEGF = c("ESR1", "PGR", "VEGFA")
Endothelial_markers = c("CD34", "SEMA3G", "GJA5",
                        "ACKR1", "PLVAP",
                        "VEGFA", "VEGFB", "ACTA2", "GUCY1A2", "IGF1", "ESR1", "PGR", "AR",
                        "PROX1", "FLT4")
Wang_Alonso_Endothelium = c("RNASE1", "PECAM1", "PCDH17", "VWF", "ADGRL4",
                            "CD34", "SEMA3G", "GJA5", "ACKR1", "PLVAP")
Schupp.Artery <- c("CLDN10", "GJA5", "GJA4", "FBLIM1", "FBLN5", "FBLN2", "MGP", "BGN", 
           "LTBP4", "LTBP1", "FN1", "SERPINE229", "CPAMD8", "CXCL12", "EFBN2", 
           "SEMA3G", "VEGFA", "NOS1", "PDE3A", "PDE4D", "DKK2", "DLL4", "HEY1", 
           "SOX5", "SOX17", "HES4", "PRDM16") 
Schupp.Capillary = c("CA416", "CYB5A", "ACVRL1/TMEM100", 
           "ADGRF5", "ADGRL2", "F2RL3", "IFNGR1", "VIPR1", "ADRB1", "ARHGAP6", 
           "IFI27", "PREX1", "PRKCE", "SGK1", "SH2D3C", "SORBS1", "PRX", "SPARC", 
           "EMP2", "ITGA1", "SLC9A3R2", "AFF3", "MEIS1") 
Schupp.Vein = c("VCAM1", "SELP", "SELE", 
           "ACKR128", "NR2F2", "ADAMTS9", "IGFBP7", "HDAC9", "RORA", "ACTN1", "LDLRAD3", 
           "LDLRAD4", "LRRC1")


# Running markerplotting on seurat clusters with RNA assay
Idents(object = endo.integrated) <- "seurat_clusters"
DefaultAssay(endo.integrated) = "RNA"
Norm.assay = "RNA"

run = Marker_plotting(marker = Alonso_Endothelial_All, name = "Alonso_Endothelial_All", group = "seurat_clusters")
run = Marker_plotting(marker = Alonso_Artery_vein, name = "Alonso_Artery_vein", group = "seurat_clusters")
run = Marker_plotting(marker = VEGF_genes, name = "VEGF_genes", group = "seurat_clusters")
run = Marker_plotting(marker = VEGF_receptor, name = "VEGF_receptor", group = "seurat_clusters")
run = Marker_plotting(marker = Hormones, name = "Hormones", group = "seurat_clusters")
run = Marker_plotting(marker = Horomones_VEGF, name = "Hormones_VEGF", group = "seurat_clusters")
run = Marker_plotting(marker = Endothelial_markers, name = "Endothelial_markers", group = "seurat_clusters")

run = Marker_plotting(marker = Wang_Alonso_Endothelium, name = "Wang_Alonso_Endothelium", group = "seurat_clusters")
run = Marker_plotting(marker = Schupp.Artery, name = "Schupp.Artery", group = "seurat_clusters")
run = Marker_plotting(marker = Schupp.Capillary, name = "Schupp.Capillary", group = "seurat_clusters")
run = Marker_plotting(marker = Schupp.Vein, name = "Schupp.Vein", group = "seurat_clusters")

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
  # Labelling all as stroma
  print("Plotting Endothelial labels")
  endo.integrated$Endothelial_labelled = endo.integrated$seurat_clusters
  endo.integrated = SetIdent(endo.integrated, value = "Endothelial_labelled")
  endo.integrated = RenameIdents(endo.integrated, "0" = "Endothelial Vein",
                                 "1" = "Lymphatic",
                                 "2" = "Endothelial Vein",
                                 "3" = "Endothelial Vein",
                                 "4" = "Endothelial Artery",
                                 "5" = "Endothelial Vein",
                                 "6" = "Endothelial proliferative",
                                 "7" = "Endothelial Vein",
                                 "8" = "Mesenchymal",
                                 "9" = "Endothelial Vein")
  endo.integrated[["Endothelial_labelled"]] = Idents(object = endo.integrated)
  DimPlot(endo.integrated, reduction = "umap", group.by = "Endothelial_labelled",
          label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir, Project_name, "_UMAP_Endothelial_labelled.pdf"), dpi = 700)
  
  # Ordering the object and adding colours
  Idents(endo.integrated) = "Endothelial_labelled"
  new.order.labels = c("Endothelial Vein", "Endothelial Artery", "Endothelial proliferative", "Mesenchymal", 
                       "Lymphatic")
  Idents(endo.integrated) <- factor(Idents(endo.integrated), levels= new.order.labels)
  endo.integrated$Endothelial_labelled <- factor(endo.integrated$Endothelial_labelled, levels= new.order.labels)
  
  # 3 colors
  Dimplot.colors = c("#8491B4CC", "#B09C8599", "#F39B7F99", "#00A08799", "#7E6148CC")
  
  # Dimplot of stroma, no groups
  Idents(endo.integrated) = "Endothelial_labelled"
  DimPlot(endo.integrated, cols = Dimplot.colors, raster = FALSE)
  ggsave2(paste0(Output.dir, Project_name, "_no-labels_UMAP_Stroma_figure_ptdefault.pdf"), dpi = 700)
  
  DimPlot(endo.integrated, cols = Dimplot.colors, raster = FALSE, pt.size = 0.5)
  ggsave2(paste0(Output.dir, Project_name, "_no-labels_UMAP_Stroma_figure_pt0.5.pdf"), dpi = 700)
  
  # Replotting the dotplots, one group
  DotPlot(object = endo.integrated, features = Endothelial_markers, group.by = "Endothelial_labelled",
          cols = c("dodgerblue", "firebrick"), dot.min = 0.1, assay = "RNA", dot.scale = 10) +
    theme(axis.text.x=element_text(angle=90, vjust = 0.3, hjust = 1))
  ggsave2(paste0(Output.dir, Project_name, "_Dotplot_Endothelial_Curated.pdf"), dpi = 700, width = 12)
  
  DotPlot(object = endo.integrated, features = Endothelial_markers, group.by = "Endothelial_labelled",
          cols = c("navajowhite", "firebrick"), dot.min = 0, dot.scale = 10) +
    theme(axis.text.x=element_text(angle=90))
  ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Endothelial_Marker_Dotplot.pdf"),
          dpi = 700, width = 14)
  
  DimPlot(endo.integrated, reduction = "umap", split.by = "Group_Stage")
  ggsave2(paste0(Output.dir, Project_name, "_UMAP_Endothelial_groups_labelled_Split-Group.pdf"), dpi = 700,
          height = 12, width = 20)
  
  # Save the labelled and re-ordered object
  SaveH5Seurat(endo.integrated, paste0(Output.dir,Project_name, "_reclustered_labelled.h5seurat"), overwrite = TRUE)
  
  ##Idents(endo.copy) <- factor(Idents(endo.copy), levels= new.order.labels)
  ###endo.copy$Immune_labelled <- factor(endo.copy$Immune_labelled, levels= new.order.labels)
  ###Dimplot.color = RColorBrewer::brewer.pal(16, "Set1")
  DimPlot(endo.integrated, cols = Dimplot.colors, raster = FALSE)
  ggsave2(paste0(Output.dir, Project_name, "_no-labels_UMAP_Stroma_grouped_figure_ptdefault.pdf"), dpi = 700)
  
  DimPlot(endo.integrated, cols = Dimplot.colors, raster = FALSE, pt.size = 0.5)
  ggsave2(paste0(Output.dir, Project_name, "_no-labels_UMAP_Stroma_grouped_figure_pt0.5.pdf"), dpi = 700)
  
  DimPlot(endo.integrated, cols = Dimplot.colors, label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir, Project_name, "_labels_UMAP_Stroma_grouped_figure.pdf"), dpi = 700)
  
  endo.average = AverageExpression(endo.integrated, return.seurat = TRUE)
  
  DoHeatmap(endo.average, features = Endothelial_markers, raster = FALSE,
            group.bar = TRUE, group.colors = c(rep("white", length(levels(endo.average)))),
            draw.lines = FALSE, angle = 45, hjust = 0) + 
    scale_fill_gradient2(low = "#2570B7", mid = "seashell", midpoint = 0, high = "#DC0000FF") +
    theme(axis.text.y = element_text(face = "italic"))
  
  ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Endothelial_groups_Marker_Heatmap.pdf"),
          dpi = 700, height = 14, width = 8)
  
  Idents(endo.average) <- factor(Idents(endo.average), levels = rev(new.order.labels))
  DoHeatmap(endo.average, features = Endothelial_markers, raster = FALSE, 
            group.bar = TRUE, group.colors = c(rep("white", length(levels(endo.average)))),
            draw.lines = FALSE, angle = 270, hjust = 1) + 
    scale_fill_gradient2(low = "#2570B7", mid = "seashell", midpoint = 0, high = "#DC0000FF") + 
    theme(axis.text.y = element_text(face = "italic", angle = 315))
  
  ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Endothelial_groups_Marker_Heatmap_flipped.pdf"), 
          dpi = 700, height = 14, width = 8)
  
  # Plotting the integrated object
  print("Plotting of UMAP")
  UMAP_Ident = DimPlot(endo.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
  ggsave2(plot = UMAP_Ident, filename = paste0(Output.dir,Project_name,"_labelled_UMAP_Ident.pdf"), dpi = 700)
  UMAP_Phase = DimPlot(endo.integrated, reduction = "umap", group.by = "Phase")
  ggsave2(plot = UMAP_Phase, filename = paste0(Output.dir,Project_name,"_UMAP_Phase.pdf"), dpi = 700)
  UMAP_Label = DimPlot(endo.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
  ggsave2(plot = UMAP_Label, filename = paste0(Output.dir,Project_name,"_labelled_UMAP_Clusters.pdf"), dpi = 700)
  UMAP_Group = DimPlot(endo.integrated, reduction = "umap", group.by = "Group_Stage",)
  ggsave2(plot = UMAP_Group, filename = paste0(Output.dir,Project_name,"_UMAP_Group_Stage.pdf"), dpi = 700)
  UMAP_original_label = DimPlot(endo.integrated, reduction = "umap", group.by = "Labelled_Clusters_SCT.1",)
  ggsave2(plot =  UMAP_original_label, filename = paste0(Output.dir,Project_name,"_Labelled_Clusters_SCT.1.pdf"), dpi = 700)
  print("Done with basic UMAP")  
  
  # QC plotting of newly labbeled clusters      
  Plot_CQ_Endothelial = Plotting_QC(x = endo.integrated, group = "Endothelial_labelled")
  
}

Plot_CQ_Group.Stage = Plotting_QC(x = endo.integrated, group = "Group_Stage")
Plot_CQ_Group.Stage = Plotting_QC(x = endo.integrated, group = "Group_Treatment")
Plot_CQ_Pat_nr = Plotting_QC(x = endo.integrated, group = "orig.ident")
Plot_CQ_seurat_clusters = Plotting_QC(x = endo.integrated, group = "seurat_clusters")

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
  Idents(object = endo.integrated) <- "seurat_clusters"
  DefaultAssay(endo.integrated) = "RNA"
  Norm.assay = "RNA"
  
  print("Finding markers")
  endo.markers <-FindAllMarkers(endo.integrated, assay = "RNA")
  print("Saving markers as .csv")
  write.csv(endo.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_seurat_clusters_FindAllMarkers.csv"), quote = F)
  print("Saving markers as .rds")
  saveRDS(endo.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_seurat_clusters_Endo_markers.rds"))
  print("Extracting top 10 markers per cluster")
  top10 = endo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  
  print("Generating heatmap")
  endo.heatmap = DoHeatmap(endo.integrated, features = top10$gene) + NoLegend()
  ggsave2(filename = paste0(Output.dir, Project_name, "_", Norm.assay, "_seurat_clusters_Top10_Genes_Heatmap.pdf"),
          plot = endo.heatmap,
          dpi = 700)
  
  # Setting cluster to be tested to SCT clustered
  Idents(object = endo.integrated) <- "Endothelial_labelled"
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
