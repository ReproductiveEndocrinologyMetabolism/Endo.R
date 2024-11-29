#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(pheatmap)

# Setting directories and flags
Project_name = "Endo_All_Stromal_uSMC"
Input.dir = "Output/4_Subsetting/"
FindMarker.flag = TRUE
ScaleData.flag = FALSE
SingleR.flag = FALSE
annotation.done = TRUE
theme_set(theme_cowplot())

# Setting the output directory
if (dir.exists(path = paste0("Output/5_Stromal_uSMC_Clustering")) == FALSE) {
  print(paste0("Generating output directory Output/5_Stromal_uSMC_Clustering"))
  dir.create(path = paste0("Output/5_Stromal_uSMC_Clustering"), recursive = TRUE)
  Output.dir = paste0("Output/5_Stromal_uSMC_Clustering/")
} else if (dir.exists(path = paste0("Output/5_Stromal_uSMC_Clustering")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/5_Stromal_uSMC_Clustering/")
} else {
  print("Error with output directory")
}

# Loading reclustered Seurat
print("Seurat object loading")
endo.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_reclustered.h5seurat"))

# Loading immune markers from FindAllMarkers
cluster.markers = readRDS("Output/4_Subsetting/Endo_All_Stromal_FindAllMarkers.rds")
cluster.markers.top10 = readRDS("Output/4_Subsetting/Endo_All_Stromal_FindAllMarkers_top10.rds")

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
Alonso_Stromal_All = c("ACTA2", "MYH11", "C7", "OGN", #FIBRO C7
                       "IGF1A", "PCOLCE", #Stroma
                       "MMP11", "CRABP2", "ECM1", #Endometrial stroma
                       "FOXO1", "IL15", "CEBPB", "PDGFRA", #Decidualized stroma
                       "AR", "ESR1", "PGR",
                       "NCAM1")

Alonso_Stromal_uSMC_All = c("ACTA2", "MYH11", "ACTG2", "RGS5", "NTRK2", "GUCY1A2", 
                       "STEAP4", "IGF1", "C7", "OGN", "MMP11", "NCAM1")
Wang_Stromal_uSMC_All = c("ACTA2", "MCAM", "BGN", "NOTH3", "GUCY1A2", "RGS5",
                          "COL5A1", "LUM", "COL6A3", "CRISPLD2", "COL6A1", "DCN1")

Alonso_fibroblast_stroma = c("ACTA2", "C7", "OGN", "IFG1", "PCOLCE", "MMP11", "CRABP2", 
                             "ECM1", "FOXO1", "IL15", "CFD", "CEBPB", "PDGFRA")
Alonso_PV = c("ACTA2", "MYH11", "ACTG2", "RGS5", "STEAP4", "GUCY1A2", "IGF1", "PCOLCE")
Wang_stromal = c("DCN", "COL6A1", "CRISPLD2", "COL6A3", "LUM", "COL5A1")
Wang_uSMC = c("ACTA2", "MCAM", "BGN", "NOTCH3", "GUCY1A2", "RGS5", "NCAM1", "CD56")
Endometriosis_Fonseca_markers = c("MME", "ESR1", "PGR", "IGF1", "MMP11", "CRABP2", "ECM1",
                                  "PAEP", "FOXO1", "IL15",
                                  "DCN", "PDGFRA", "THY1", "COL1A1", "PCOLCE",
                                  "C7", "OGN", "CFD",
                                  "ITGB1", "FAP", "ACTA2",
                                  "ACTG2", "MCAM",
                                  "GAS5")
Hormones = c("AR", "ESR1", "PGR")

uSMC_PV_marker = c("ACTA2", "MYH11", "ACTG2", "GUCY1A2", "NCAM1", "MCAM", "BGN", "NOTCH3", "RGS5")
Stroma_markers = c("")

# Running markerplotting on seurat clusters with RNA assay
Idents(object = endo.integrated) <- "seurat_clusters"
DefaultAssay(endo.integrated) = "RNA"
Norm.assay = "RNA"

run = Marker_plotting(marker = Alonso_Stromal_uSMC_All, name = "Alonso_Stromal_uSMC_All", group = "seurat_clusters")
run = Marker_plotting(marker = Wang_Stromal_uSMC_All, name = "Wang_Stromal_uSMC_All", group = "seurat_clusters")
run = Marker_plotting(marker = Alonso_Stromal_All, name = "Alonso_Stromal_All", group = "seurat_clusters")
run = Marker_plotting(marker = Alonso_fibroblast_stroma, name = "Alonso_fibroblast_stroma", group = "seurat_clusters")
run = Marker_plotting(marker = Alonso_PV, name = "Alonso_PV", group = "seurat_clusters")
run = Marker_plotting(marker = Wang_stromal, name = "Wang_stromal", group = "seurat_clusters")
run = Marker_plotting(marker = Wang_uSMC, name = "Wang_uSMC", group = "seurat_clusters")
run = Marker_plotting(marker = Hormones, name = "Hormones", group = "seurat_clusters")
run = Marker_plotting(marker = Endometriosis_Fonseca_markers, name = "Endometriosis_Fonseca_markers", group = "seurat_clusters")

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
  print("Plotting Stromal labels")
  endo.integrated$Stromal_labelled = endo.integrated$seurat_clusters
  endo.integrated = SetIdent(endo.integrated, value = "Stromal_labelled")
  endo.integrated = RenameIdents(endo.integrated, "0" = "Stroma 1",
                                 "1" = "Stroma 1",
                                 "2" = "Stroma 1",
                                 "3" = "Stroma 1",
                                 "4" = "Stroma 2",
                                 "5" = "Stroma 1",
                                 "6" = "Stroma 1",
                                 "7" = "Fibroblast",
                                 "8" = "Stroma 1",
                                 "9" = "Stroma 1",
                                 "10" = "Stroma 1",
                                 "11" = "Stroma 1",
                                 "12" = "Fibroblast",
                                 "13" = "Stroma 1",
                                 "14" = "Stroma 1",
                                 "15" = "Stroma 1",
                                 "16" = "uSMC",
                                 "17" = "Stroma proliferative",
                                 "18" = "Stroma 1",
                                 "19" = "Stroma 1",
                                 "20" = "Stroma 1",
                                 "21" = "Stroma 1",
                                 "22" = "Stroma 1",
                                 "23" = "Stroma proliferative",
                                 "24" = "Stroma 1",
                                 "25" = "Stroma 2")
  endo.integrated[["Stromal_labelled"]] = Idents(object = endo.integrated)
  DimPlot(endo.integrated, reduction = "umap", group.by = "Stromal_labelled",
          label = TRUE, repel = TRUE)
  ggsave2(paste0(Output.dir, Project_name, "_UMAP_Stromal_labelled.pdf"), dpi = 700)
  
  # Ordering the object and adding colours
  Idents(endo.integrated) = "Stromal_labelled"
  new.order.labels = c("Stroma 1", "Stroma 2", "Stroma proliferative", "Fibroblast", "uSMC")
  Idents(endo.integrated) <- factor(Idents(endo.integrated), levels= new.order.labels)
  endo.integrated$Stromal_labelled <- factor(endo.integrated$Stromal_labelled, levels= new.order.labels)
  
  # 3 colors
  #Dimplot.colors = c("#E4211C", "#FAA0A1", "#F29403", "#91D1C2CC") #Cellchat
  Dimplot.colors = c("#E64B35CC", "#DC000099", "#FAA0A1", "#91D1C2CC", "#F39B7FCC") #NPG
  
  # Dimplot of stroma, no groups
  Idents(endo.integrated) = "Stromal_labelled"
  DimPlot(endo.integrated, cols = Dimplot.colors, raster = FALSE)
  ggsave2(paste0(Output.dir, Project_name, "_no-labels_UMAP_Stroma_figure_ptdefault.pdf"), dpi = 700)
  
  DimPlot(endo.integrated, cols = Dimplot.colors, raster = FALSE, pt.size = 0.5)
  ggsave2(paste0(Output.dir, Project_name, "_no-labels_UMAP_Stroma_figure_pt0.5.pdf"), dpi = 700)
  
  # Replotting the dotplots, one group
  DotPlot(object = endo.integrated, features = Alonso_Stromal_All, group.by = "Stromal_labelled",
          cols = c("dodgerblue", "firebrick"), dot.min = 0.1, assay = "RNA", dot.scale = 10) +
    theme(axis.text.x=element_text(angle=90, vjust = 0.3, hjust = 1))
  ggsave2(paste0(Output.dir, Project_name, "_Dotplot_Stromal_Curated.pdf"), dpi = 700, width = 12)
  
  DotPlot(object = endo.integrated, features = Alonso_Stromal_All, group.by = "Stromal_labelled",
          cols = c("navajowhite", "firebrick"), dot.min = 0, dot.scale = 10) +
    theme(axis.text.x=element_text(angle=90))
  ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Stromal_Marker_Dotplot.pdf"),
          dpi = 700, width = 14)
  
  DimPlot(endo.integrated, reduction = "umap", split.by = "Group_Stage")
  ggsave2(paste0(Output.dir, Project_name, "_UMAP_Stromal_groups_labelled_Split-Group.pdf"), dpi = 700,
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
  
  DoHeatmap(endo.average, features = Endometriosis_Fonseca_markers, raster = FALSE,
            group.bar = TRUE, group.colors = c(rep("white", length(levels(endo.average)))),
            draw.lines = FALSE, angle = 45, hjust = 0) + 
    scale_fill_gradient2(low = "#2570B7", mid = "seashell", midpoint = 0, high = "#DC0000FF") +
    theme(axis.text.y = element_text(face = "italic"))

  ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Stromal_groups_Marker_Heatmap.pdf"),
          dpi = 700, height = 14, width = 8)
  
  Idents(endo.average) <- factor(Idents(endo.average), levels = rev(new.order.labels))
  DoHeatmap(endo.average, features = Endometriosis_Fonseca_markers, raster = FALSE, 
            group.bar = TRUE, group.colors = c(rep("white", length(levels(endo.average)))),
            draw.lines = FALSE, angle = 270, hjust = 1) + 
    scale_fill_gradient2(low = "#2570B7", mid = "seashell", midpoint = 0, high = "#DC0000FF") + 
    theme(axis.text.y = element_text(face = "italic", angle = 315))
  
  ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Stromal_groups_Marker_Heatmap_flipped.pdf"), 
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
  Plot_CQ_Stromal = Plotting_QC(x = endo.integrated, group = "Stromal_labelled")

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
  Idents(object = endo.integrated) <- "Stromal_labelled"
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