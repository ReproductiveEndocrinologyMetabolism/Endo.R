#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(pheatmap)

# Set the variables
Project.name = "Endo_10x_snRNA_All_QC"
#seurat.x = "Output/Plotting_Average_Heatmap/Endo_All_Combined_labels_Control_PCOS.h5seurat"
seurat.x = "Output/9_Main_relabelled/Endo_All_Combined_labels.h5seurat"
Group.cols = c("#A0A0A0", "#D098B2", "#95BFE1", "#65D46E")
Group.order = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS")
Sample.cols = rev(c(rep("#A0A0A0", 5), rep("#D098B2", 12), rep("#95BFE1", 7), rep("#65D46E", 3)))
Sample.order = rev(c("1_Control", "2_Control", "3_Control", "4_Control", "5_Control",
                     "6_PCOS_Lifestyle_W0", "7_PCOS_Lifestyle_W0", "8_PCOS_Lifestyle_W0", "9_PCOS_EA_W0",
                     "10_PCOS_Metformin_W0", "11_PCOS_Metformin_W0",  "11_PCOS_Metformin_W0",
                     "12_PCOS_Metformin_W0",  "13_PCOS_Metformin_W0", "14_PCOS_Metformin_W0",
                     "15_PCOS_Metformin_W0", "16_PCOS_Metformin_W0",  "17_PCOS_Metformin_W0",
                     "18_PCOS_Lifestyle_W16", "19_PCOS_Lifestyle_W16", "20_PCOS_Lifestyle_W16",
                     "21_PCOS_Metformin_W16", "22_PCOS_Metformin_W16", "23_PCOS_Metformin_W16",
                     "24_PCOS_Metformin_W16", "25_PCOS_Metformin_W16", "26_PCOS_Metformin_W16",
                     "27_PCOS_Metformin_W16")) # Lifestyle

main.markers.col = c(rep("#4DBBD5CC", 3), rep("#E64B35CC", 3), rep("#F39B7FCC", 3),
                     rep("#91D1C2CC", 3), rep("#00A087CC", 4), rep("#8491B4CC", 2),
                     rep("#7E6148CC", 2))
main.markers = c("EPCAM", "CPM", "LGR5", # Epithelium
                    "IGF1", "DCN", "COL6A1", # Stroma
                    "GUCY1A2", "ACTA2", "NOTCH3", # uSMC
                    "CD14", "CSF1R", "LYZ",  #Myeloids
                    "STK17B", "NCAM1", "CCL5", "CD2", # Lymphoids
                    "PCDH17", "VWF", # Endothelial
                    "PROX1", "FLT4") # Lymphatic
names(main.markers.col) = main.markers

Main.CT.order = c("Stromal", "Epithelium", "uSMC", "Myeloid", "Lymphoid", "Endothelial", "Lymphatic", "Undefined")
Main.CT.cols = list("#E64B35CC", # Stromal RED
                  "#4DBBD5CC", # Epithelium BLUE
                  "#F39B7FCC", # uSMC RED-PINKISH
                  "#91D1C2CC", # Myeloid
                  "#00A087CC", # Lymphoid GREEN
                  "#8491B4CC", # Endothelial Grey
                  "#7E6148CC", # Lymphatic Brown
                  "black") # Undefined
names(Main.CT.cols) = Main.CT.order

subtype.markers = c("PTGS1", "VTCN1", "SLC26A7", # Lumenal
                    "LGR5", "KRT5", "WNT7A", # Sox9+ LGR5+
                    "CPM", "IHH", "EMID1", "PPARG", # Sox9+ LGR5-
                    "C2CD4A", "SLC18A2", "CXCL14", # SOX9 prolif
                    "MKI67", "HMGB2","AR", # AR+
                    "CDC20B", "CCNO", "HES6", # Ciliated
                    "ESR1", "PGR", "IGF1", "ECM1", "PAEP", #Stroma
                    "OGN", "TOP2A", # Stroma proliferative
                    "THY1", "COL1A1", "PCOLCE", "C7", # Fibroblast
                    "ACTA2", "ACTG2", "MCAM", # uSMC
                    "FOXP3", "CD3G", "CD27", "IKZF4", "CD8A", "CCL5", "CD4", 
                    "IL7R", "RORC", #Treg, T-cell CD8, T-cell CD4, ILC3
                    "MS4A1", "IGHM", #B-celler
                    "NCAM1", "ITGA1", "SPINK2", "CSF1", "CD160", "GNLY", # uNK 1-3
                    "IL3RA", "LILRA4", "PLD4",# pDC
                    "EBI3", "CCR7", "CCL19", # Migratory DC
                    "BATF3", "CADM1", "CLEC9A", #DC1
                    "CLEC10A", "FCER1A", "CD1C", #DC2
                    "CD14", "SELENOP", "HMOX1", "IL1B", #uM 1-2
                    "CPA3", "KIT", "MS4A2", # Mast cells
                    "PECAM1", "CD34" , # Endothelial
                    "ACKR1", "PLVAP", # Endothelial Vein
                    "SEMA3G", "GJA5", # Endothelial Artery
                    "COL3A1", "WNT5A", "MMP11", # Mesenchymal
                    "PROX1", "FLT4") # Lymphatic

main_subtype.markers = c("EPCAM", # Epithelium main
                         "PTGS1", "VTCN1", "SLC26A7", # Lumenal
                    "LGR5", "KRT5", "WNT7A", # Sox9+ LGR5+
                    "CPM", "IHH", "EMID1", "PPARG", # Sox9+ LGR5-
                    "C2CD4A", "SLC18A2", "CXCL14", # SOX9 prolif
                    "MKI67", "HMGB2","AR", # AR+
                    "CDC20B", "CCNO", "HES6", # Ciliated,
                    "IGF1", "DCN", "COL6A1", # Stroma main
                    "ESR1", "PGR", "ECM1", "PAEP", #Stroma
                    "OGN", "TOP2A", # Stroma proliferative
                    "THY1", "COL1A1", "PCOLCE", "C7", # Fibroblast
                    "GUCY1A2", "NOTCH3", "ACTA2", "ACTG2", "MCAM", # uSMC
                    "STK17B", "NCAM1", "CCL5", "CD2", # Lymphoid main
                    "FOXP3", "CD3G", "CD27", "IKZF4", "CD8A", "CD4", 
                    "IL7R", "RORC", #Treg, T-cell CD8, T-cell CD4, ILC3
                    "MS4A1", "IGHM", #B-celler
                    "ITGA1", "SPINK2", "CSF1", "CD160", "GNLY", # uNK 1-3
                    "CD14", "CSF1R", "LYZ",  #Myeloids main
                    "IL3RA", "LILRA4", "PLD4",# pDC
                    "EBI3", "CCR7", "CCL19", # Migratory DC
                    "BATF3", "CADM1", "CLEC9A", #DC1
                    "CLEC10A", "FCER1A", "CD1C", #DC2
                    "SELENOP", "HMOX1", "IL1B", #uM 1-2
                    "CPA3", "KIT", "MS4A2", # Mast cells
                    "PCDH17", "VWF", # Endothelial
                    "PECAM1", "CD34" , # Endothelial
                    "ACKR1", "PLVAP", # Endothelial Vein
                    "SEMA3G", "GJA5", # Endothelial Artery
                    "WNT5A", "MMP11", # Mesenchymal
                    "PROX1", "FLT4") # Lymphatic


subtype.order = c("Lumenal", "SOX9+ LGR5+", "SOX9+ LGR5-",              
                  "SOX9+ proliferative", "AR+", "Ciliated", # Epithelium
                  "Stroma 1", "Stroma 2", "Stroma proliferative", 
                  "Fibroblast", "uSMC", # Stroma
                  "Tregs", "T-cells CD8+", "T-cells CD4+",
                  "ILC3", "B-cells", "uNK 1", "uNK 2", "uNK 3",
                  "pDC", "Migratory DC", "DC1", "DC2", 
                  "uM 1", "uM 2", "Mast cells", # Immune
                  "Endothelial Vein", "Endothelial Artery", "Endothelial proliferative", 
                  "Mesenchymal", "Lymphatic", NA) # Endothelial


# Set the output directory
if (dir.exists(path = "Output/Paper_Plotting_QC") == FALSE) {
  print("Output/Paper_Plotting_QC")
  dir.create(path = "Output/Paper_Plotting_QC", recursive = TRUE)
  Output.dir = "Output/Paper_Plotting_QC/"
} else if (dir.exists(path = "Output/Paper_Plotting_QC") == TRUE) {
  print("Directory exists")
  Output.dir = "Output/Paper_Plotting_QC/"
} else {
  print("Error with output directory")
}

### Plot out quality control values of the main Seurat object

# Load the seurat object
seurat.x = LoadH5Seurat(seurat.x)

# Set assay to RNA
DefaultAssay(seurat.x) = "RNA"

# Order idents
seurat.x$Labelled_Clusters_fine = factor(seurat.x$Labelled_Clusters_fine, levels = Main.CT.order)
seurat.x$orig.ident = factor(seurat.x$orig.ident, levels = Sample.order)
seurat.x$Group_Stage = factor(seurat.x$Group_Stage, levels = Group.order)
seurat.x$Combined_labels = factor(seurat.x$Combined_labels, levels = subtype.order)

DimPlot(object = seurat.x, cols = Main.CT.cols, group.by = "Labelled_Clusters_fine") + theme_map()
ggsave(filename=paste0(Output.dir, Project.name,"_Dimplot_Cluster_with_Undefined.pdf"), dpi = 300)

DimPlot(seurat.x, group.by = "orig.ident") + theme_map()
ggsave(filename=paste0(Output.dir, Project.name,"_Dimplot_SampleID.pdf"), dpi = 300)

DimPlot(seurat.x, group.by = "Group_Stage", cols = Group.cols) + theme_map()
ggsave(filename=paste0(Output.dir, Project.name,"_Dimplot_Group_Stage.pdf"), dpi = 300)

DimPlot(seurat.x, group.by = "seurat_clusters") + theme_map()
ggsave(filename=paste0(Output.dir, Project.name,"_Dimplot_seurat_clusters.pdf"), dpi = 300)

DimPlot(seurat.x, group.by = "Phase") + theme_map()
ggsave(filename=paste0(Output.dir, Project.name,"_Dimplot_Phase.pdf"), dpi = 300)

# Plot out the main markers
for (marker.x in main.markers) {
  
  marker.col = main.markers.col[marker.x]

  FeaturePlot(object = seurat.x, features = marker.x, cols = c("#FFFFCC", marker.col), reduction = "umap", max.cutoff = 3) + theme_map()
  ggsave(filename=paste0(Output.dir, Project.name, "_", marker.x, "_Featureplot_Beige_Max3.pdf"), dpi = 300)
  
  FeaturePlot(object = seurat.x, features = marker.x, cols = c("#FFFFCC", "firebrick"), reduction = "umap", max.cutoff = 3) + theme_map()
  ggsave(filename=paste0(Output.dir, Project.name, "_", marker.x, "_Featureplot_Default_Max3.pdf"), dpi = 300)
  
}

# Subset Undefined before QC
Idents(object = seurat.x) <- "Combined_labels"
seurat.x = subset(seurat.x, idents = subtype.order[-32])
seurat.x$Combined_labels = factor(seurat.x$Combined_labels, levels = rev(subtype.order))
seurat.x = NormalizeData(seurat.x)


seurat.x$Combined_labels = factor(seurat.x$Combined_labels, levels = rev(subtype.order))

DotPlot(object = seurat.x, features = main_subtype.markers,
        cols = c("#FFFFCC", "firebrick"), group.by = "Combined_labels", 
        scale = FALSE, dot.min = 0.1, cluster.idents = FALSE) + 
  theme(axis.text.x=element_text(angle=45, hjust=1,  face = "italic"))
ggsave(filename=paste0(Output.dir, Project.name, "_DotPlot_Combined_labels_Subtype_Markers_Log2.pdf"), 
       width = 21, height = 10, dpi = 300)

DotPlot(object = seurat.x, features = main_subtype.markers,
        cols = c("#FFFFCC", "firebrick"), group.by = "Combined_labels", 
        scale = TRUE, dot.min = 0.1, cluster.idents = FALSE) + 
  theme(axis.text.x=element_text(angle=45, hjust=1,  face = "italic"))
ggsave(filename=paste0(Output.dir, Project.name, "_DotPlot_Combined_labels_Subtype_Markers_Scaled.pdf"), 
       width = 21, height = 10, dpi = 300)

seurat.x$Combined_labels = factor(seurat.x$Combined_labels, levels = subtype.order)

DotPlot(object = seurat.x, features = rev(main_subtype.markers),
        cols = c("#FFFFCC", "firebrick"), group.by = "Combined_labels", 
        scale = FALSE, dot.min = 0.1, cluster.idents = FALSE) + theme(axis.text.x=element_text(angle=45, hjust=1), 
                               axis.text.y=element_text(angle=0, hjust=1, face = "italic")) + coord_flip()
ggsave(filename=paste0(Output.dir, Project.name, "_DotPlot_Combined_labels_Subtype_Markers_Log2_flip.pdf"), 
       width = 10, height = 18, dpi = 300)

DotPlot(object = seurat.x, features = rev(main_subtype.markers),
        cols = c("#FFFFCC", "firebrick"), group.by = "Combined_labels", 
        scale = TRUE, dot.min = 0.1, cluster.idents = FALSE) + theme(axis.text.x=element_text(angle=45, hjust=1), 
                                                                      axis.text.y=element_text(angle=0, hjust=1, face = "italic")) + coord_flip()
ggsave(filename=paste0(Output.dir, Project.name, "_DotPlot_Combined_labels_Subtype_Markers_Scaled_flip.pdf"), 
       width = 10, height = 18, dpi = 300)

# Generate a averaged seurat
seurat.average = LoadH5Seurat(file = "Output/Paper_Plotting_QC/Endo_10x_snRNA_All_QC_Combined_Labels_Averaged.h5seurat")
seurat.average$orig.ident = factor(seurat.x$orig.ident, levels = rev(subtype.order))

# Do a heatmap of the markers
average.df = seurat.average@assays$RNA$scale.data[rownames(seurat.average@assays$RNA$scale.data) %in% main_subtype.markers, ]
average.df = average.df[match(main_subtype.markers, rownames(average.df)), ]

plot.x = pheatmap::pheatmap(average.df, cluster_rows = FALSE, cluster_cols = FALSE,
                            cellwidth = 10, cellheight = 5, fontsize = 5)
ggsave(plot = plot.x, filename=paste0(Output.dir, Project.name, "_Heatmap_Combined_labels_Subtype_Markers_Scaled.pdf"), 
       width = 8, height = 8, dpi = 300)

plot.x = pheatmap::pheatmap(average.df, cluster_rows = FALSE, cluster_cols = FALSE,
                            cellwidth = 5, cellheight = 5, fontsize = 5, gaps_col = c(6, 11, 19, 26))
ggsave(plot = plot.x, filename=paste0(Output.dir, Project.name, "_Heatmap_Combined_labels_Subtype_Markers_Scaled_Gaps_col.pdf"), 
       width = 8, height = 8, dpi = 300)

plot.x = pheatmap::pheatmap(average.df, cluster_rows = FALSE, cluster_cols = TRUE,
                            cellwidth = 10, cellheight = 5, fontsize = 5)
ggsave(plot = plot.x, filename=paste0(Output.dir, Project.name, "_Heatmap_Combined_labels_Subtype_Markers_Scaled_Clustered.pdf"), 
       width = 8, height = 8, dpi = 300)

# Do plotting of QC values
VlnPlot(seurat.x, features = "nFeature_RNA", cols = rev(Sample.cols), group.by = "orig.ident", pt.size = 0) + theme_cowplot() +
  theme(legend.position = 'none')
ggsave(filename=paste0(Output.dir, Project.name, "_VlnPlot_nFeature_SampleID.pdf"), dpi = 300)

VlnPlot(seurat.x, features = "nCount_RNA", cols = rev(Sample.cols), group.by = "orig.ident", pt.size = 0) + theme_cowplot() +
  theme(legend.position = 'none')
ggsave(filename=paste0(Output.dir, Project.name, "_VlnPlot_nCount_SampleID.pdf"), dpi = 300)

seurat.x$Group_Stage = factor(seurat.x$Group_Stage, levels = rev(Group.order))
RidgePlot(seurat.x, features = "nFeature_RNA", cols = rev(Group.cols), group.by = "Group_Stage") + theme_cowplot() +
  theme(legend.position = 'none')
ggsave(filename=paste0(Output.dir, Project.name, "_Rigdeplot_nFeature_Group.pdf"), dpi = 300)

RidgePlot(seurat.x, features = "nCount_RNA", cols = rev(Group.cols), group.by = "Group_Stage") + theme_cowplot() +
  theme(legend.position = 'none')
ggsave(filename=paste0(Output.dir, Project.name, "_Rigdeplot_nCount_Group.pdf"), dpi = 300)

# Plot Ridgeplots of QC values
seurat.x$orig.ident = factor(seurat.x$orig.ident, levels = Sample.order)
print(levels(seurat.x$orig.ident))

RidgePlot(seurat.x, features = "nFeature_RNA", cols = Sample.cols, group.by = "orig.ident") + theme_cowplot() +
  theme(legend.position = 'none')
ggsave(filename=paste0(Output.dir, Project.name, "_Rigdeplot_nFeature_SampleID.pdf"), dpi = 300)

RidgePlot(seurat.x, features = "nCount_RNA", cols = Sample.cols, group.by = "orig.ident") + theme_cowplot() +
  theme(legend.position = 'none')
ggsave(filename=paste0(Output.dir, Project.name, "_Rigdeplot_nCount_SampleID.pdf"), dpi = 300)

seurat.x$Group_Stage = factor(seurat.x$Group_Stage, levels = rev(Group.order))
RidgePlot(seurat.x, features = "nFeature_RNA", cols = rev(Group.cols), group.by = "Group_Stage") + theme_cowplot() +
  theme(legend.position = 'none')
ggsave(filename=paste0(Output.dir, Project.name, "_Rigdeplot_nFeature_Group.pdf"), dpi = 300)

RidgePlot(seurat.x, features = "nCount_RNA", cols = rev(Group.cols), group.by = "Group_Stage") + theme_cowplot() +
  theme(legend.position = 'none')
ggsave(filename=paste0(Output.dir, Project.name, "_Rigdeplot_nCount_Group.pdf"), dpi = 300)


