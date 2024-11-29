#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(future)

#setwd("/mnt/data/guseri/10x_analysi/All_sampless")
Project_name = "Endo_All"
FindMarker.flag = TRUE
Norm.assay = "RNA"
Input.dir = "Output/1_Integrated/"
ScaleData.flag = TRUE
Ident.group = "seurat_clusters"

if (dir.exists(path = paste0("Output/2_Clustering")) == FALSE) {
  print(paste0("Generating output directory Output/2_Clustering"))
  dir.create(path = paste0("Output/2_Clustering"), recursive = TRUE)
  Output.dir = paste0("Output/2_Clustering/")
} else if (dir.exists(path = paste0("Output/2_Clustering")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/2_Clustering/")
} else {
  print("Error with output directory")
}

# Plotting function
Marker_plotting <- function(marker, name, group = Ident.group, selected.assay = Norm.assay) {
  
  print(paste("Plotting", name, "with", selected.assay))
  
  DefaultAssay(endo.integrated) = selected.assay
  VlnPlot(endo.integrated, features = marker, pt.size = 0, sort = "increasing", group.by = group)
  ggsave2(paste0(Output.dir, Project_name, "_", group, "_", name, "_", Norm.assay, "_Marker_ViolinPlot.pdf"), dpi = 700)
  
  FeaturePlot(endo.integrated, features = marker)
  ggsave2(paste0(Output.dir, Project_name, "_", group, "_", name, "_", Norm.assay, "_Marker_FeaturePlot.pdf"), dpi = 700)
  
  DotPlot(object = endo.integrated, features = marker, group.by = group) + RotatedAxis()
  ggsave2(paste0(Output.dir, Project_name, "_", group, "_", name, "_", Norm.assay, "_Marker_Dotplot.pdf"), dpi = 700)
}

# Marker list
# Updated cell cluster markers
Stromal_cells = c("IGF1", "GUCY1A2", "AR", "MMP11", "ECM1", "CFD", "CEBPB")
Epithelium = c("SCGB2A2", "CPM", "ESR1", "HMGB2", "AR", "PGR", "EPCAM", "KRT8", "LGR5")
Lumenal_epithelium = c("PGR", "ESR1", "CPM", "HMGB2", "LGR5")
Glandular_epithelium = c("SCGB2A2", "CPM", "ESR1", "HMGB2", "AR")
Ciliated_epithelium = c("MMP7", "HMGB2", "MUC12", "CDC20B", "CCNO", "HES6")
Glandular_secretory_epithelium = c("PAEP", "C2CD4A", "SLC18A2", "CXCL14", "GPX3")
Smooth_Muscle = c("ACTA2", "MYH11", "ACTG2", "NCAM1", "GUCY1A2")
Perivascular = c("ACTA2", "MYH11", "NTRK2", "GUCY1A2") # Cluster 17
Endothelial_Artery_Vein = c("CD34", "GJA5", "SEMA3G") # Cluster 15
Lymphatic = c("PROX1", "FLT4") # Cluster 8, 15, 17
Lymphoid = c("PTPRC", "NCAM1", "CD3G", "ACTA2") # Cluster 11, 14
B.cells = c("CD38", "IGF1") #Cluster 11, 14
Myeloid = c("PTPRC", "CD14", "CSF1R") # Cluster 15
Epithelial = c("KRT8", "KRT18",	"EPCAM","CD24")
Non.decidualised_endometrial.stromal.cells = c("ACTA2", "MYH11", "GUCY1A2", "IGF1")
monocytes = c("CD14", "CD11b", "CCR2","CD16")
mast_cells = c("KIT", "MS4A2")
stroma_2 = c("DCN", "COL6A3", "LUM")
epithelium_2 = c("EPCAM", "WFDC2", "KLF5", "UCA1", "TACSTD2")
ciliated_2 = c("FOXJ1", "C9orf24", "CDHR3", "DYDC2")
lymphocytes_2 = c("PTPRC", "CCL5", "STK17B")
macrophages_2 = c("LYZ", "HLA-DQA1", "AIF1", "MS4A6A")
endothelial_2 = c("RNASE1", "PECAM1", "PCDH17", "VWF", "ADGRL4")
SMC_2 = c("ACTA2", "MCAM", "NOTCH3", "GUCY1A2", "RGS5")
Perivascular_2 = c("STEAP4", "MYH1")

## Use the top 3 markers
stroma_3 = c("DCN", "COL6A3", "LUM")
epithelium_3 = c("EPCAM", "WFDC2", "KLF5", "TACSTD2")
ciliated_3 = c("FOXJ1", "C9orf24", "CDHR3")
lymphocytes_3 = c("PTPRC", "CCL5", "STK17B")
macrophages_3 = c("LYZ", "HLA-DQA1", "MS4A6A")
endothelial_3 = c("PCDH17", "VWF", "ADGRL4")
SMC_3 = c("ACTA2", "NOTCH3", "GUCY1A2")
Perivascular_3 = c("STEAP4", "MYH1")

# Loading Seurat
print("Seurat object loading")
endo.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_", Norm.assay, "_labeled_integrated.h5seurat"))

# Setting cluster to be tested to SCT clustered
DefaultAssay(endo.integrated) = "RNA"
Idents(object = endo.integrated) <- "seurat_clusters"

#Plotting SCT clustering UMAP
print("Plotting of UMAP")
UMAP_Ident = DimPlot(endo.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_labelled_UMAP_Ident.pdf"), dpi = 700)
UMAP_Phase = DimPlot(endo.integrated, reduction = "umap", group.by = "Phase")
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_UMAP_Phase.pdf"), dpi = 700)
UMAP_Label = DimPlot(endo.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_labelled_UMAP_Clusters.pdf"), dpi = 700)
UMAP_Label = DimPlot(endo.integrated, reduction = "umap", group.by = "seurat_clusters", label = FALSE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_UMAP_Clusters.pdf"), dpi = 700)
UMAP_Group = DimPlot(endo.integrated, reduction = "umap", group.by = "Group_Stage")
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_UMAP_Group_Stage.pdf"), dpi = 700)
print("Done with basic UMAP of log2 normalised object")

# Running marker plotting
run = Marker_plotting(marker = Stromal_cells, name = "Stromal_cells")
run = Marker_plotting(marker = Epithelium, name = "Epithelium")
run = Marker_plotting(marker = Lumenal_epithelium, name = "Lumenal_epithelium")
run = Marker_plotting(marker = Glandular_epithelium, name = "Glandular_epithelium")
run = Marker_plotting(marker = Ciliated_epithelium, name = "Ciliated_epithelium")
run = Marker_plotting(marker = Glandular_secretory_epithelium, name = "Glandular_secretory_epithelium")
run = Marker_plotting(marker = Smooth_Muscle, name = "Smooth_Muscle")
run = Marker_plotting(marker = Perivascular, name = "Perivascular")
run = Marker_plotting(marker = Endothelial_Artery_Vein, name = "Endothelial_Artery_Vein")
run = Marker_plotting(marker = Lymphatic, name = "Lymphatic")
run = Marker_plotting(marker = Lymphoid, name = "Lymphoid")
run = Marker_plotting(marker = B.cells, name = "B.cells")
run = Marker_plotting(marker = Myeloid, name = "Myeloid")
run = Marker_plotting(marker = Epithelial, name = "Epithelial")
run = Marker_plotting(marker = Non.decidualised_endometrial.stromal.cells, name = "Non.decidualised_endometrial.stromal.cells")
run = Marker_plotting(marker = monocytes, name = "monocytes")
run = Marker_plotting(marker = mast_cells, name = "mast_cells")
run = Marker_plotting(marker = stroma_2, name = "stroma_2")
run = Marker_plotting(marker = epithelium_2, name = "epithelium_2")
run = Marker_plotting(marker = ciliated_2, name = "ciliated_2")
run = Marker_plotting(marker = lymphocytes_2, name = "lymphocytes_2")
run = Marker_plotting(marker = macrophages_2, name = "macrophages_2")
run = Marker_plotting(marker = endothelial_2, name = "endothelial_2")
run = Marker_plotting(marker = SMC_2, name = "SMC_2")
run = Marker_plotting(marker = Perivascular_2, name = "Perivascular_2")
run = Marker_plotting(marker = stroma_3, name = "stroma_3")
run = Marker_plotting(marker = epithelium_3, name = "epithelium_3")
run = Marker_plotting(marker = ciliated_3, name = "ciliated_3")
run = Marker_plotting(marker = lymphocytes_3, name = "lymphocytes_3")
run = Marker_plotting(marker = macrophages_3, name = "macrophages_3")
run = Marker_plotting(marker = endothelial_3, name = "endothelial_3")
run = Marker_plotting(marker = SMC_3, name = "SMC_3")
run = Marker_plotting(marker = Perivascular_3, name = "Perivascular_3")

# Run FindAllMarkers to identify markers for clusters
if (FindMarker.flag == TRUE) {
  
  # Data is re-scaled after subsetting as the mean and SD will have changed 
  if (ScaleData.flag == TRUE) {
    Norm.assay = "RNA"
    print("Setting DefaultAssay to RNA and log normalising it")
    DefaultAssay(endo.integrated) = "RNA"
    endo.integrated = NormalizeData(endo.integrated)
    endo.integrated = FindVariableFeatures(endo.integrated)
    all.genes = rownames(endo.integrated)
    endo.integrated = ScaleData(endo.integrated, features = all.genes)  #endo.integrated = FindVariableFeatures(endo.integrated)
    
  } else if (ScaleData.flag == FALSE) {
    print("Setting DefaultAssay to SCT")
    DefaultAssay(endo.integrated) = "SCT" #Alternative is RNA
  }
  
  SaveH5Seurat(endo.integrated, paste0(Input.dir, Project_name, "_", Norm.assay, "_labeled_integrated.h5seurat"), overwrite = TRUE)

  print("Finding markers")
  endo.markers <-FindAllMarkers(endo.integrated, assay = Norm.assay, 
                                logfc.threshold = 0.5, min.pct = 0.25)
  print("Saving markers as .csv")
  write.csv(endo.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_FindAllMarkers.csv"), quote = F)
  print("Saving markers as .rds")
  saveRDS(endo.markers, paste0(Output.dir, Project_name, "_", Norm.assay, "_Endo_markers.rds"))
  print("Extracting top 10 markers per cluster")
  top10 = endo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) # Maybe change to p-value?
  
  
  print("Generating heatmap")
  endo.heatmap = DoHeatmap(endo.integrated, features = top10$gene) + NoLegend()
  ggsave2(filename = paste0(Output.dir, Project_name, "_", Norm.assay, "_Top10_Genes_Heatmap.pdf"),
          plot = endo.heatmap,
          dpi = 700)
}

if (SingleR.flag == TRUE) {
  
  library(celldex)
  library(SingleR)
  
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
  ggsave2(paste0(Output.dir,Project_name,"_labelled_SingleR_Monaco_main.pdf"), dpi = 700)
  
  # Testing and plotting DbImmune main
  endo.SingleR$immune_labels_DbImmune_main = singleR.main.DbImmune$pruned.labels
  Idents(endo.SingleR) = "immune_labels_DbImmune_main"
  DimPlot(endo.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_SingleR_DbImmune_main.pdf"), dpi = 700)
  
  # Testing and plotting Blueprint main
  endo.SingleR$immune_labels_Blueprint_main = singleR.main.Blueprint$pruned.labels
  Idents(endo.SingleR) = "immune_labels_Blueprint_main"
  DimPlot(endo.SingleR, label = T , repel = T)
  ggsave2(paste0(Output.dir,Project_name,"_labelled_SingleR_Blueprint_main.pdf"), dpi = 700)
}

print("Done")

