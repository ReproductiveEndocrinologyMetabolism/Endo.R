#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(future)

# Setting the parameters
Project_name = "Endo_All"
FindMarker.flag = FALSE
Input.dir = "Output/1_Integrated/"
SCT.flag = TRUE
Log2.flag = FALSE

if (dir.exists(path = paste0("Output/3_Labeling")) == FALSE) {
  print(paste0("Generating output directory Output/3_Labeling"))
  dir.create(path = paste0("Output/3_Labeling"), recursive = TRUE)
  Output.dir = paste0("Output/3_Labeling/")
} else if (dir.exists(path = paste0("Output/3_Labeling")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/3_Labeling/")
} else {
  print("Error with output directory")
}

# Plotting QC metrics 
Plotting_QC <- function(x, group = "") {
  
  # Generating ridgeplots
  RidgePlot(endo.integrated, features = "percent.mt", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_mtDNA.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "percent.ribo", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_ribo.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "percent.hb", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_hb.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "nCount_RNA", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_nCount_RNA.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "nFeature_RNA", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_nFeature_RNA.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "S.Score", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_S-Score.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "G2M.Score", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_G2M-Score.pdf"), dpi = 700)
  
  # Generating violin plots
  VlnPlot(endo.integrated, features = "percent.mt", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_mtDNA.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "percent.ribo", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_ribo.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "percent.hb", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_hb.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "nCount_RNA", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_nCount_RNA.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "nFeature_RNA", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_nFeature_RNA.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "S.Score", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_S-Score.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "G2M.Score", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_G2M-Score.pdf"), dpi = 700)
  
  # Make stacked barplot of group population per cell type to determine proportunality
  # between Control and PCOS in each cell type
  
  Idents(object = endo.integrated) <- group
  pt <- table(Idents(endo.integrated), endo.integrated$Group_Stage)
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)
  
  # Add colours for the colour blind
  barplot_prop = ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Group") +
    ylab("Proportion") +
    theme(legend.title = element_blank()) +
    theme_cowplot()
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_barplot_proportions_Stage.pdf"), dpi = 700)
  
  #Make stacked barplot of group population per cell type to determine proportunality
  # between Control and PCOS in each cell type
  
  Idents(object = endo.integrated) <- group
  pt <- table(Idents(endo.integrated), endo.integrated$Group_Treatment)
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)
  
  # Add colours for the colour blind
  barplot_prop = ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Group") +
    ylab("Proportion") +
    theme(legend.title = element_blank()) +
    theme_cowplot()
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_barplot_proportions_Treatment.pdf"), dpi = 700)
  
}

# Loading SCT clustered Seurat object
Norm.assay = "RNA"
print("SCT clusters seurat object loading")
endo.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_", Norm.assay, "_labeled_integrated.h5seurat"))

print("Setting DefaultAssay to SCT")
# Set the default assay to RNA
DefaultAssay(endo.integrated) = "RNA" 

# Setting cluster to be tested to SCT clustered
Idents(object = endo.integrated) <- "seurat_clusters"

##### ALL SAMPLES #####
# Renaming for groups and stage
print("Plotting 1")
endo.integrated$Labelled_Clusters_SCT.1 = endo.integrated$seurat_clusters
endo.integrated = SetIdent(endo.integrated, value = "Labelled_Clusters_SCT.1")
endo.integrated = RenameIdents(endo.integrated, "0" = "Stromal",
                               "1" = "Stromal",
                               "2" = "Stromal",
                               "3" = "Epithelium",
                               "4" = "Epithelium",
                               "5" = "Epithelium",
                               "6" = "Stromal",
                               "7" = "Epithelium",
                               "8" = "Stromal",
                               "9" = "Epithelium",
                               "10" = "Epithelium",
                               "11" = "Epithelium",
                               "12" = "Epithelium",
                               "13" = "Immune",
                               "14" = "Immune",
                               "15" = "Stromal",
                               "16" = "Epithelium",
                               "17" = "Immune",
                               "18" = "Stromal",
                               "19" = "Stromal",
                               "20" = "uSMC",
                               "21" = "Epithelium",
                               "22" = "Epithelium",
                               "23" = "Endothelial",
                               "24" = "Immune",
                               "25" = "Undefined",
                               "26" = "Epithelium",
                               "27" = "Lymphatic",
                               "28" = "Immune",
                               "29" = "Epithelium")
endo.integrated[["Labelled_Clusters_SCT.1"]] = Idents(object = endo.integrated)
DimPlot(endo.integrated, reduction = "umap", group.by = "Labelled_Clusters_SCT.1", 
        label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir, Project_name, "_", Norm.assay, "_UMAP_Labelled_Clusters_1.pdf"), dpi = 700)

# Labelling on smaller resolution
print("Plotting 2")
endo.integrated$Labelled_Clusters_SCT.2 = endo.integrated$seurat_clusters
endo.integrated = SetIdent(endo.integrated, value = "Labelled_Clusters_SCT.2")
endo.integrated = RenameIdents(endo.integrated, "0" = "Stromal",
                               "1" = "Stromal",
                               "2" = "Stromal",
                               "3" = "Epithelium",
                               "4" = "Epithelium",
                               "5" = "Epithelium",
                               "6" = "Stromal",
                               "7" = "Epithelium",
                               "8" = "Stromal",
                               "9" = "Epithelium",
                               "10" = "Epithelium",
                               "11" = "Epithelium",
                               "12" = "Epithelium",
                               "13" = "Immune_Myeloid",
                               "14" = "Immune_Lymphoid_NK",
                               "15" = "Stromal",
                               "16" = "Epithelium",
                               "17" = "Immune_Lymphoid_T.cells",
                               "18" = "Stromal",
                               "19" = "Stromal",
                               "20" = "uSMC",
                               "21" = "Epithelium",
                               "22" = "Epithelium",
                               "23" = "Endothelial",
                               "24" = "Immune_cycling_lymphocytes",
                               "25" = "Immune_ILC",
                               "26" = "Epithelium",
                               "27" = "Lymphatic",
                               "28" = "Immune_mast.cells",
                               "29" = "Epithelium")
endo.integrated[["Labelled_Clusters_SCT.2"]] = Idents(object = endo.integrated)
DimPlot(endo.integrated, reduction = "umap", group.by = "Labelled_Clusters_SCT.2",
        label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir, Project_name, "_", Norm.assay, "_UMAP_Labelled_Clusters_2.pdf"), dpi = 700)

# Clean labelling of fine labels, specifically immune
# Labelling on smaller resolution
print("Plotting 3")
endo.integrated$Labelled_Clusters_fine = endo.integrated$seurat_clusters
endo.integrated = SetIdent(endo.integrated, value = "Labelled_Clusters_fine")
endo.integrated = RenameIdents(endo.integrated, "0" = "Stromal",
                               "1" = "Stromal",
                               "2" = "Stromal",
                               "3" = "Epithelium",
                               "4" = "Epithelium",
                               "5" = "Epithelium",
                               "6" = "Stromal",
                               "7" = "Epithelium",
                               "8" = "Stromal",
                               "9" = "Epithelium",
                               "10" = "Epithelium",
                               "11" = "Epithelium",
                               "12" = "Epithelium",
                               "13" = "Myeloid",
                               "14" = "Lymphoid",
                               "15" = "Stromal",
                               "16" = "Epithelium",
                               "17" = "Lymphoid",
                               "18" = "Stromal",
                               "19" = "Stromal",
                               "20" = "uSMC",
                               "21" = "Epithelium",
                               "22" = "Epithelium",
                               "23" = "Endothelial",
                               "24" = "Lymphoid",
                               "25" = "Undefined",
                               "26" = "Epithelium",
                               "27" = "Lymphatic",
                               "28" = "Myeloid",
                               "29" = "Epithelium")
endo.integrated[["Labelled_Clusters_fine"]] = Idents(object = endo.integrated)
DimPlot(endo.integrated, reduction = "umap", group.by = "Labelled_Clusters_fine",
        label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir, Project_name, "_", Norm.assay, "_UMAP_Labelled_Clusters_fine.pdf"), dpi = 700)

# QC plotting
Plot_CQ_labeled_clusters = Plotting_QC(x = endo.integrated, group = "Labelled_Clusters_fine")
Plot_CQ_labeled_clusters = Plotting_QC(x = endo.integrated, group = "Labelled_Clusters_SCT.1")
Plot_CQ_labeled_clusters = Plotting_QC(x = endo.integrated, group = "seurat_clusters")

# Log2 normalising
Idents(object = endo.integrated) <- "Labelled_Clusters_SCT.1"
DefaultAssay(endo.integrated) = "RNA"
endo.integrated = NormalizeData(endo.integrated)
endo.integrated = FindVariableFeatures(endo.integrated)
all.genes = rownames(endo.integrated)
endo.integrated = ScaleData(endo.integrated, features = all.genes)

# Save labeled object
SaveH5Seurat(object = endo.integrated, filename = paste0(Output.dir,Project_name, "_", Norm.assay,"_celltypes_integrated.h5seurat"), overwrite = TRUE)

if (FindMarker.flag == TRUE) {
  
  # Setting cluster to be tested to SCT clustered
  Idents(object = endo.integrated) <- "seurat_clusters"
  
  print("Finding markers")
  endo.markers <-FindAllMarkers(endo.integrated, assay = "RNA", logfc.threshold = 0.5, min.pct = 0.5)
  print("Saving markers as .csv")
  write.csv(endo.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_FindAllMarkers.csv"), quote = F)
  print("Saving markers as .rds")
  saveRDS(endo.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_Endo_markers.rds"))
  print("Extracting top 10 markers per cluster")
  top5 = endo.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  
  
  print("Generating heatmap")
  endo.heatmap = DoHeatmap(endo.integrated, features = top5$gene) + NoLegend()
  ggsave2(filename = paste0(Output.dir, Project_name, "_", Norm.assay, "_Top5_Genes_Heatmap.pdf"),
          plot = endo.heatmap,
          dpi = 700)
  
  # Setting cluster to be tested to SCT clustered
  Idents(object = endo.integrated) <- "Labelled_Clusters_SCT.1"
  
  print("Finding markers")
  endo.markers <-FindAllMarkers(endo.integrated, assay = "RNA", logfc.threshold = 0.5, min.pct = 0.5)
  print("Saving markers as .csv")
  write.csv(endo.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_FindAllMarkers_labelled.csv"), quote = F)
  print("Saving markers as .rds")
  saveRDS(endo.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_Endo_markers_labelled.rds"))
  print("Extracting top 10 markers per cluster")
  top5 = endo.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  
  
  print("Generating heatmap")
  endo.heatmap = DoHeatmap(endo.integrated, features = top5$gene) + NoLegend()
  ggsave2(filename = paste0(Output.dir, Project_name, "_", Norm.assay, "_Top5_Genes_Heatmap_labelled.pdf"),
          plot = endo.heatmap,
          dpi = 700)
  
}