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
library(dplyr)

# Set the variables
Input.dir = "Output/1_Label_Transfer_Annotation/" # "Output/1_Label_Transfer_Annotation/" "Output/1_Manual_Sample_Annotation/"
Seurat.label = "predicted.celltype" #"Labelled_cluster" "predicted.celltype"
Seurat.pattern = ".*UMAP_integrated_labelled_filtered.h5seurat" # "UMAP_integrated_labelled_filtered.h5seurat ".*bin30_labelled.h5Seurat"

# Set celltype order and colours
celltype.order = c("Epithelium", "Stromal", "uSMC", "Immune", "Endothelial") # c("Epithelium", "Stromal", "uSMC", "Immune", "Endothelial")
celltype.color = list("#4DBBD5CC", # Epithelium BLUE
                      "#E64B35CC", # Stromal RED
                      "#F39B7FCC", # uSMC RED-PINKISH
                      "#00A087CC", # Immune GREEN
                      "#8491B4CC") # Endothelial Grey
names(celltype.color) = celltype.order

# Load the seurat object
seurat.list = list.files(Input.dir, pattern = Seurat.pattern, full.names = TRUE, recursive = TRUE)
print(seurat.list)

Plotting.main <- function(seurat.x, ct.order = celltype.order, ct.col = celltype.color, 
                          seurat.ident = Seurat.label, set.max = 3, automated.labeling = TRUE) {
  
  # Get the sample name
  seurat.sample = gsub(".*Output/.*//(.*?)/.*", "\\1", seurat.x)
  print(seurat.sample)
  
  if (automated.labeling == FALSE) {
    # Generate output directory
    if (dir.exists(path = paste0("Output/4_Main_plotting_Manual/", seurat.sample)) == FALSE) {
      print(paste0("Output/4_Main_plotting_Manual/", seurat.sample))
      dir.create(path =paste0("Output/4_Main_plotting_Manual/", seurat.sample), recursive = TRUE)
      Output.dir = paste0("Output/4_Main_plotting_Manual/", seurat.sample, "/")
    } else if (dir.exists(path = paste0("Output/4_Main_plotting_Manual/", seurat.sample)) == TRUE) {
      print("Directory exists")
      Output.dir = paste0("Output/4_Main_plotting_Manual/", seurat.sample, "/")
    } else {
      print("Error with output directory")
    }
  } else if (automated.labeling == TRUE) {
    # Generate output directory
    if (dir.exists(path = paste0("Output/4_Main_plotting_Automated/", seurat.sample)) == FALSE) {
      print(paste0("Output/4_Main_plotting_Automated/", seurat.sample))
      dir.create(path =paste0("Output/4_Main_plotting_Automated/", seurat.sample), recursive = TRUE)
      Output.dir = paste0("Output/4_Main_plotting_Automated/", seurat.sample, "/")
    } else if (dir.exists(path = paste0("Output/4_Main_plotting_Automated/", seurat.sample)) == TRUE) {
      print("Directory exists")
      Output.dir = paste0("Output/4_Main_plotting_Automated/", seurat.sample, "/")
    } else {
      print("Error with output directory")
    }
  }
  
  # Load the seurat object
  seurat.x = LoadH5Seurat(seurat.x)
  
  # Set the ident with selected label
  Idents(seurat.x) = seurat.ident
  
  # Reorder seurat based on selected order
  Idents(seurat.x) <- factor(Idents(seurat.x), levels= ct.order)
  
  # Save the selected ident for future plotting
  seurat.rds = seurat.x[[seurat.ident]]
  saveRDS(seurat.rds, file = paste0(Output.dir,seurat.sample, "_Labelled_Cluster.rds"))
  seurat.rds = NULL
  
  levels(seurat.x)
  ct.col = ct.col[names(ct.col) %in% levels(seurat.x)]
  
  # Do a dimplot with the clusters
  DimPlot(seurat.x, cols = ct.col, reduction = "spatial", pt.size = 0.1) + theme_cowplot()
  ggsave2(paste0(Output.dir,seurat.sample, "_Dimplot_main_clusters.pdf"))
  
  DimPlot(seurat.x, cols = ct.col, reduction = "spatial", pt.size = 0.1) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Dimplot_main_clusters_mini.pdf"))
  
  # Do a dimplot with the clusters
  DimPlot(seurat.x, cols = ct.col, reduction = "umap", pt.size = 0.1) + theme_cowplot()
  ggsave2(paste0(Output.dir,seurat.sample, "_Dimplot_main_clusters_UMAP.pdf"))
  
  DimPlot(seurat.x, cols = ct.col, reduction = "umap", pt.size = 0.1) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Dimplot_main_clusters_mini_UMAP.pdf"))
  
  # Do a proportion barplot and return the data
  # Check and plot proportions of celltypes. Save the output
  prop.df = data.frame(prop.table(table(seurat.x[[seurat.ident]]))*100)
  colnames(prop.df) = c("Celltype", "Percentage")
  prop.df$Sample = seurat.sample

  # Reorder the dataframe based on celltype order
  prop.df$Celltype <- factor(prop.df$Celltype, levels = ct.order)
  prop.df <- prop.df[order(prop.df$Celltype), ]
  
  # Save the dataframe
  saveRDS(prop.df, file = paste0(Output.dir, seurat.sample, "_celltype_proportions.rds"))
  
  # Plot additional QC variables
  FeaturePlot(seurat.x, features = "nFeature_RNA", reduction = "spatial", cols = c("#FFFFCC", "#B10026"), pt.size = 0.1) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_nFeature.pdf"), dpi = 300)
  
  if (automated.labeling == TRUE) {
    FeaturePlot(seurat.x, features = "predicted.celltype.score", reduction = "spatial", cols = c("#FFFFCC", "#0099FF"), pt.size = 0.1)
    ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_prediction_score.pdf"), dpi = 300)
  }
  
  ##### Plot key markers with celltype specific colours
  FeaturePlot(seurat.x, features = c("EPCAM", "CPM", "LGR5"), reduction = "spatial", cols = c("#FFFFCC", "#4DBBD5CC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_epithelium_markers.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "EPCAM", reduction = "spatial", cols = c("#FFFFCC", "#4DBBD5CC"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_EPCAM.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "EPCAM", reduction = "umap", cols = c("#FFFFCC", "#4DBBD5CC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_UMAP_EPCAM.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("IGF1", "DCN", "COL6A1"), reduction = "spatial", cols = c("#FFFFCC", "#E64B35CC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_stromal_markers.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "IGF1", reduction = "umap", cols = c("#FFFFCC", "#E64B35CC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_UMAP_IGF1.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("GUCY1A2", "ACTA2", "NOTCH3"), reduction = "spatial", cols = c("#FFFFB2", "#F39B7FCC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_uSMC_markers.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "ACTA2", reduction = "spatial", cols = c("#FFFFCC", "#F39B7FCC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_ACTA2.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("CD14", "CSF1R", "LYZ"), reduction = "spatial", cols = c("#FFFFB2", "#00A087CC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_lymphoid_markers.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "CD14", reduction = "spatial", cols = c("#FFFFCC", "#00A087CC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_CD14.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("NCAM1", "CCL5", "CD2"), reduction = "spatial", cols = c("#FFFFCC", "#91D1C2CC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_myeloid_markers.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "CD2", reduction = "spatial", cols = c("#FFFFCC", "#91D1C2CC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_CD2.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("PCDH17", "VWF", "PROX1", "FLT4"), reduction = "spatial", cols = c("#FFFFCC", "#8491B4CC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_endothelial_markers.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "VWF", reduction = "spatial", cols = c("#FFFFCC", "#8491B4CC"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_VWF.pdf"), dpi = 300)
  
  ##### Plot key markers with consistent colours
  FeaturePlot(seurat.x, features = c("EPCAM", "CPM", "LGR5"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_epithelium_markers_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "EPCAM", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_EPCAM_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("IGF1", "DCN", "COL6A1"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_stromal_markers_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "IGF1", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_IGF1_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("GUCY1A2", "ACTA2", "NOTCH3"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_uSMC_markers_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "ACTA2", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_ACTA2_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("CD14", "CSF1R", "LYZ"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_lymphoid_markers_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "CD14", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_CD14_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("NCAM1", "CCL5", "CD2"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_myeloid_markers_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "CD2", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_CD2_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("PCDH17", "VWF", "PROX1", "FLT4"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_endothelial_markers_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "VWF", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_VWF_monocolor.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "VWF", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_VWF_monocolor.pdf"), dpi = 300)
  
  
  VlnPlot(seurat.x, features = c("EPCAM", "IGF1", "ACTA2", "CD14", "CD2", "VWF"), pt.size = 0) + theme_cowplot()
  ggsave2(paste0(Output.dir,seurat.sample, "_Vlnplot_keymarkers.pdf"), dpi = 300)
  
  # Minimal theme
  ##### Plot key markers with consistent colours
  FeaturePlot(seurat.x, features = c("EPCAM", "CPM", "LGR5"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max)  + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_epithelium_markers_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "EPCAM", reduction = "spatial", cols = c("#FFFFCC", "#4DBBD5CC"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_EPCAM_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("IGF1", "DCN", "COL6A1"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_stromal_markers_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "IGF1", reduction = "spatial", cols = c("#FFFFCC", "#E64B35CC"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_IGF1_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("GUCY1A2", "ACTA2", "NOTCH3"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_uSMC_markers_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "ACTA2", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_ACTA2_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("CD14", "CSF1R", "LYZ"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_lymphoid_markers_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "CD14", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_CD14_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("NCAM1", "CCL5", "CD2"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_myeloid_markers_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "CD2", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_CD2_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = c("PCDH17", "VWF", "PROX1", "FLT4"), reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_endothelial_markers_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "VWF", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_VWF_monocolor_mini.pdf"), dpi = 300)
  
  FeaturePlot(seurat.x, features = "VWF", reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_Featureplot_VWF_monocolor_mini.pdf"), dpi = 300)
  
  
  VlnPlot(seurat.x, features = c("EPCAM", "IGF1", "ACTA2", "CD14", "CD2", "VWF"), pt.size = 0) + theme_cowplot()
  ggsave2(paste0(Output.dir,seurat.sample, "_Vlnplot_keymarkers.pdf"), dpi = 300)
  
  
  # Return the prop.df for downstream plotting
  return(prop.df)

}

# Run the main plotting function
main.output = lapply(seurat.list, Plotting.main)
saveRDS(main.output, file = "Output/4_Main_plotting/Proportion_tables_automated_labeling.rds")

# Plot out celltype proportions of multple groups
main.output = readRDS(file = "Output/4_Main_plotting/Proportion_tables_automated_labeling.rds")
prop.df.all = bind_rows(main.output[2:7])
samples.df = unique(prop.df.all$Sample)
sample.order = c("02_Ctrl_bin30", "28_Ctrl_bin30", "08_PCOS_bin30", "11_PCOS_bin30", "20_LS_bin30", "22_Met_bin30")
prop.df.all$Sample <- factor(prop.df.all$Sample, levels = sample.order)
prop.df.all <- prop.df.all[order(prop.df.all$Sample), ]

ggplot(prop.df.all, aes(x = Sample, y = Percentage, fill = Celltype)) +
  geom_bar(stat = "identity") + theme_cowplot() +
  theme(legend.position = "right") +
  labs(x = "Sample", y = "Percentage", fill = "Cell Type") +
  scale_fill_manual(values = celltype.color)
ggsave2(filename = "Output/4_Main_plotting/Proportions_barplot_automated.pdf", 
        dpi = 300)

# Function to plot specific region
Plot_region <- function(seurat.x, cord.x, cord.y, marker.plot, output.name, ct.sample,  
                        ct.order = celltype.order, ct.col = celltype.color,
                        seurat.ident = Seurat.label, set.pt.size = 2, set.max = 5,
                        do.subset = FALSE) {
  
  # Get the sample name
  seurat.sample = gsub(".*Output/.*//(.*?)/.*", "\\1", seurat.x)
  print(seurat.sample)
  
  # Generate output directory
  if (dir.exists(path = paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name)) == FALSE) {
    print(paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name))
    dir.create(path =paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name), recursive = TRUE)
    Output.dir = paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name, "/")
  } else if (dir.exists(path = paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name)) == TRUE) {
    print("Directory exists")
    Output.dir = paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name, "/")
  } else {
    print("Error with output directory")
  }
  
  # Load the seurat object
  seurat.x = LoadH5Seurat(seurat.x)
  
  # Set the ident with selected label
  Idents(seurat.x) = seurat.ident
  
  # Reorder seurat based on selected order
  Idents(seurat.x) <- factor(Idents(seurat.x), levels= ct.order)
  
  levels(seurat.x)
  ct.col = ct.col[names(ct.col) %in% levels(seurat.x)]
  
  # Plot the marker on the full spatial object
  for (marker.x in names(marker.plot)) {

    FeaturePlot(seurat.x, features = marker.x, reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + 
      theme_map() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.title = element_blank())
    ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_", marker.x, "_Featureplot_full_spatial_monocolor_mini.pdf"), dpi = 300)
    
    FeaturePlot(seurat.x, features = marker.x, reduction = "spatial", cols = c("#FFFFCC", marker.plot[[marker.x]]), pt.size = 0.1, max.cutoff = set.max) + 
      theme_map() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.title = element_blank())
    ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_", marker.x, "_Featureplot_full_spatial_ctColor_mini.pdf"), dpi = 300)
    
  }
  
  
  # Plot the gene markers on only subset labelled regions if TRUE (Default)
  if (do.subset == TRUE) {
    
    seurat.ct = subset(seurat.x, idents = ct.sample)
    DimPlot(seurat.ct, cols = ct.col[[ct.sample]], reduction = "spatial", pt.size = 1) + theme_map()
    ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_", ct.sample, "_Dimplot_subset_labelled.pdf"))
    
    for (marker.x in names(marker.plot)) {
    
      FeaturePlot(seurat.ct, features = marker.x, reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = 0.1, max.cutoff = set.max) + 
        theme_map() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_blank())
      ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_", ct.sample, "_", marker.x, "_Featureplot_subset_labelled_monocolor_mini.pdf"), dpi = 300)
      
      FeaturePlot(seurat.ct, features = marker.x, reduction = "spatial", cols = c("#FFFFCC", marker.plot[[marker.x]]), pt.size = 0.1, max.cutoff = set.max) + 
        theme_map() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_blank())
      ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_", ct.sample, "_", marker.x, "_Featureplot_subset_labelled_ctColor_mini.pdf"), dpi = 300)
      
      
    }
    
    # Remove the object to save space
    seurat.ct = NULL
    
  }
  
  # Illustrade subset area
  # Do a dimplot with the clusters
  DimPlot(seurat.x, cols = ct.col, reduction = "spatial", pt.size = 1) + 
    geom_hline(yintercept=c(cord.y), linetype='dashed', col = 'black') +
    geom_vline(xintercept=c(cord.x), linetype='dashed', col = 'black') +
    theme_cowplot()
  ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_Dimplot_main_clusters_dashed_region.pdf"))
  
  # Subset the spatial dimension
  seurat.x = subset(seurat.x, subset = spatial_1 > cord.x[1] & spatial_1 < cord.x[2] & 
                            spatial_2 > cord.y[1] & spatial_2 < cord.y[2])
  
  # Do a dimplot with the clusters
  DimPlot(seurat.x, cols = ct.col, reduction = "spatial", pt.size = set.pt.size) + theme_map()
  ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_Dimplot_main_clusters_selected_region_mini.pdf"))
  
  # Plot out selected genes
  for (marker.x in names(marker.plot)) {
    
    FeaturePlot(seurat.x, features = marker.x, reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = set.pt.size, max.cutoff = set.max) + 
      theme_map() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.title = element_blank())
    ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_", marker.x, "_Featureplot_selected_region_monocolor_mini.pdf"), dpi = 300)
    
    FeaturePlot(seurat.x, features = marker.x, reduction = "spatial", cols = c("#FFFFCC", marker.plot[[marker.x]]), pt.size = set.pt.size, max.cutoff = set.max) + 
      theme_map() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.title = element_blank())
    ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_", marker.x, "_Featureplot_selected_region_ctColor_mini.pdf"), dpi = 300)
    
  }
  
  # If subset is True, subset the selected region for the target cell type
  if (do.subset == TRUE) {
    seurat.x = subset(seurat.x, idents = ct.sample)
    DimPlot(seurat.x, cols = ct.col[[ct.sample]], reduction = "spatial", pt.size = set.pt.size) + theme_map()
    ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_Dimplot_subset_labelled_selected_region.pdf"))
    
    marker.plot.list = list()
    for (marker.x in names(marker.plot)) {
      # Plot the markers on the subsetted selected region
      FeaturePlot(seurat.x, features = marker.x, reduction = "spatial", cols = c("#FFFFCC", "#CC0033"), pt.size = set.pt.size, max.cutoff = set.max) + 
        theme_map() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_blank())
      ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_", ct.sample, "_", marker.x, "_Featureplot_subset_selected_region_monocolor_mini.pdf"), dpi = 300)
      
      plot.x = FeaturePlot(seurat.x, features = marker.x, reduction = "spatial", cols = c("#FFFFCC", marker.plot[[marker.x]]), pt.size = set.pt.size, max.cutoff = set.max) + 
        theme_map() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_blank())
      ggsave2(plot = plot.x, filename = paste0(Output.dir,seurat.sample, "_", output.name, "_", ct.sample, "_", marker.x, "_Featureplot_subset_selected_region_ctColor_mini.pdf"), dpi = 300)
      
      marker.plot.list[[marker.x]] = plot.x
      
    }
    
    return(marker.plot.list) 
    
  }
  
}

# Plotting epithelium region
# Plotting all markers, with the same colors
epithelium.markers = list(EPCAM = "#4DBBD5CC", ESR1 = "#4DBBD5CC", # Epithelium
                          PTGS1 = "#65C2A5", VTCN1 = "#65C2A5", SLC26A7 = "#65C2A5", LGR5 = "#65C2A5", #Lumenal
                          KRT5 = "#D3020D", WNT7A = "#D3020D", #SOX9+ LGR5+
                          CPM = "#2570B7", IHH = "#2570B7", EMID1 = "#2570B7", PPARG = "#2570B7", # SOX9+ LGR5-
                          C2CD4A = "#FAA0A1", SLC18A2 = "#FAA0A1", PAEP = "#FAA0A1", CXCL14 = "#FAA0A1", 
                          MKI67 = "#FAA0A1", HMGB2 = "#FAA0A1", # SOX9+ proliferative
                          AR = "#FF7F00", # AR+
                          CDC20B = "#DD8D61", CCNO = "#DD8D61", HES6 = "#DD8D61") # Ciliated #DD8D61

#epithelium.markers = c("EPCAM", 'ESR1', 'LGR5', 'WNT7A', 'IHH', 'AR')
sample_08_PCOS = Plot_region(seurat.x = seurat.list[2], cord.x = c(17000, 19000), cord.y = c(13000, 15000), marker.plot = epithelium.markers, output.name = "Epithelium", ct.sample = "Epithelium")
sample_20_LS = Plot_region(seurat.x = seurat.list[3], cord.x = c(9000, 11000), cord.y = c(11000, 13000), marker.plot = epithelium.markers, output.name = "Epithelium", ct.sample = "Epithelium")
sample_11_PCOS = Plot_region(seurat.x = seurat.list[4], cord.x = c(14000, 16000), cord.y = c(12000, 14000), marker.plot = epithelium.markers, output.name = "Epithelium", ct.sample = "Epithelium")
sample_22_Met = Plot_region(seurat.x = seurat.list[5], cord.x = c(15000, 17000), cord.y = c(14000, 16000), marker.plot = epithelium.markers, output.name = "Epithelium", ct.sample = "Epithelium")
sample_02_Ctrl = Plot_region(seurat.x = seurat.list[6], cord.x = c(11500, 13500), cord.y = c(9500, 11500), marker.plot = epithelium.markers, output.name = "Epithelium", ct.sample = "Epithelium")
sample_28_Ctrl = Plot_region(seurat.x = seurat.list[7], cord.x = c(12000, 14000), cord.y = c(8000, 10000), marker.plot = epithelium.markers, output.name = "Epithelium", ct.sample = "Epithelium")


# Plotting stroma region
stromal.markers = c(ESR1 = "#E64B35CC", PGR = "#E64B35CC", IGF1 = "#E64B35CC", ECM1 = "#E64B35CC", # Stroma 1
                    PAEP = "#4DBBD599", # Stroma 2
                    OGN = "#FAA0A1", TOP2A = "#FAA0A1", MKI67 = "#FAA0A1", # Stroma proliferative
                    THY1 = "#91D1C2CC", COL1A1 = "#91D1C2CC", PCOLCE = "#91D1C2CC", C7 = "#91D1C2CC", # Fibroblast
                    ACTA2 = "#F39B7FCC", ACTG2 = "#F39B7FCC", MCAM = "#F39B7FCC", # uSMC
                    EPCAM = "#4DBBD5CC", AR = "#E64B35CC") # Sanity check

sample_08_PCOS = Plot_region(seurat.x = seurat.list[2], cord.x = c(17000,19000), cord.y = c(14000, 16000), marker.plot = stromal.markers, output.name = "Stromal", ct.sample = "Stromal")
sample_20_LS = Plot_region(seurat.x = seurat.list[3], cord.x = c(11000, 13000), cord.y = c(10000, 12000), marker.plot = stromal.markers, output.name = "Stromal", ct.sample = "Stromal")
sample_11_PCOS = Plot_region(seurat.x = seurat.list[4], cord.x = c(15000, 17000), cord.y = c(11000, 13000), marker.plot = stromal.markers, output.name = "Stromal", ct.sample = "Stromal")
sample_22_Met = Plot_region(seurat.x = seurat.list[5], cord.x = c(13000, 15000), cord.y = c(13000, 15000), marker.plot = stromal.markers, output.name = "Stromal", ct.sample = "Stromal")
sample_02_Ctrl = Plot_region(seurat.x = seurat.list[6], cord.x = c(15000, 17000), cord.y = c(10000, 12000), marker.plot = stromal.markers, output.name = "Stromal", ct.sample = "Stromal")
sample_28_Ctrl = Plot_region(seurat.x = seurat.list[7], cord.x = c(14000, 16000), cord.y = c(10000, 12000), marker.plot = stromal.markers, output.name = "Stromal", ct.sample = "Stromal")

# Plotting DEG and Hormones
DEG_Hormone.markers = list(ESR1 = "#CC0033", PGR = "#CC0033", AR = "#CC0033", # Hormones
                           ROBO2 = "#CC0033", ITGA2 = "#CC0033", LAMC1 = "#CC0033", 
                           ADAMTS9 = "#CC0033", CD44 = "#CC0033", # MMP3 = "#CC0033",
                           HIF1A = "#CC0033", RUNX1 = "#CC0033", COL1A2 = "#CC0033") # DEG's

sample_08_PCOS = Plot_region(seurat.x = seurat.list[2], cord.x = c(17000, 19000), cord.y = c(13000, 15000), marker.plot = DEG_Hormone.markers, output.name = "DEG_Hormone_Max5", ct.sample = "Epithelium")
sample_20_LS = Plot_region(seurat.x = seurat.list[3], cord.x = c(9000, 11000), cord.y = c(11000, 13000), marker.plot = DEG_Hormone.markers, output.name = "DEG_Hormone_Max5", ct.sample = "Epithelium")
sample_11_PCOS = Plot_region(seurat.x = seurat.list[4], cord.x = c(14000, 16000), cord.y = c(12000, 14000), marker.plot = DEG_Hormone.markers, output.name = "DEG_Hormone_Max5", ct.sample = "Epithelium")
sample_22_Met = Plot_region(seurat.x = seurat.list[5], cord.x = c(15000, 17000), cord.y = c(14000, 16000), marker.plot = DEG_Hormone.markers, output.name = "DEG_Hormone_Max5", ct.sample = "Epithelium")
sample_02_Ctrl = Plot_region(seurat.x = seurat.list[6], cord.x = c(11500, 13500), cord.y = c(9500, 11500), marker.plot = DEG_Hormone.markers, output.name = "DEG_Hormone_Max5", ct.sample = "Epithelium")
sample_28_Ctrl = Plot_region(seurat.x = seurat.list[7], cord.x = c(12000, 14000), cord.y = c(8000, 10000), marker.plot = DEG_Hormone.markers, output.name = "DEG_Hormone_Max5", ct.sample = "Epithelium")

# Plotting DEG and Hormones
Main_10x.markers = list(EPCAM = "#4DBBD5CC", 
                        IGF1 = "#E64B35CC", 
                        GUCY1A2 = "#F39B7FCC",
                        CD14 = "#91D1C2CC", 
                        NCAM1 = "#00A087CC", 
                        VWF = "#8491B4CC", 
                        PROX1 = "#7E6148CC")

sample_08_PCOS = Plot_region(seurat.x = seurat.list[2], cord.x = c(17000, 19000), cord.y = c(13000, 15000), marker.plot = Main_10x.markers, output.name = "Main_10x_markers_Max2", ct.sample = "Epithelium", set.max = 2)
sample_20_LS = Plot_region(seurat.x = seurat.list[3], cord.x = c(9000, 11000), cord.y = c(11000, 13000), marker.plot = Main_10x.markers, output.name = "Main_10x_markers_Max2", ct.sample = "Epithelium", set.max = 2)
sample_11_PCOS = Plot_region(seurat.x = seurat.list[4], cord.x = c(14000, 16000), cord.y = c(12000, 14000), marker.plot = Main_10x.markers, output.name = "Main_10x_markers_Max2", ct.sample = "Epithelium", set.max = 2)
sample_22_Met = Plot_region(seurat.x = seurat.list[5], cord.x = c(15000, 17000), cord.y = c(14000, 16000), marker.plot = Main_10x.markers, output.name = "Main_10x_markers_Max2", ct.sample = "Epithelium", set.max = 2)
sample_02_Ctrl = Plot_region(seurat.x = seurat.list[6], cord.x = c(11500, 13500), cord.y = c(9500, 11500), marker.plot = Main_10x.markers, output.name = "Main_10x_markers_Max2", ct.sample = "Epithelium", set.max = 2)
sample_28_Ctrl = Plot_region(seurat.x = seurat.list[7], cord.x = c(12000, 14000), cord.y = c(8000, 10000), marker.plot = Main_10x.markers, output.name = "Main_10x_markers_Max2", ct.sample = "Epithelium", set.max = 2)


# Plotting CellChat ligands and receptors
# Ligand markers
CellChat_ligand_epithelium.markers = c(COL1A1 = "#4DBBD5CC", COL1A2 = "#4DBBD5CC", COL4A1 = "#4DBBD5CC", COL4A2 = "#4DBBD5CC", COL6A3 = "#4DBBD5CC", FN1 = "#4DBBD5CC", LAMA2 = "#4DBBD5CC", LAMC1 = "#4DBBD5CC", LAMC2 = "#4DBBD5CC", LAMC3 = "#4DBBD5CC",
                            SPP1 = "#4DBBD5CC", GDF7 = "#4DBBD5CC", IGF1 = "#4DBBD5CC", SEMA3E = "#4DBBD5CC", SLIT2 = "#4DBBD5CC",
                            CADM1 = "#4DBBD5CC", CNTN1 = "#4DBBD5CC", NRXN3 = "#4DBBD5CC")

CellChat_ligand_stroma.markers = c(COL1A1 = "#E64B35CC", COL1A2 = "#E64B35CC", COL4A1 = "#E64B35CC", COL4A2 = "#E64B35CC", COL6A3 = "#E64B35CC", FN1 = "#E64B35CC", LAMA2 = "#E64B35CC", LAMC1 = "#E64B35CC", LAMC2 = "#E64B35CC", LAMC3 = "#E64B35CC",
                                   SPP1 = "#E64B35CC", GDF7 = "#E64B35CC", IGF1 = "#E64B35CC", SEMA3E = "#E64B35CC", SLIT2 = "#E64B35CC",
                                   CADM1 = "#E64B35CC", CNTN1 = "#E64B35CC", NRXN3 = "#E64B35CC")

# Receptor markers
CellChat_receptor_epithelium.markers = c(ITGA1 = "#4DBBD5CC", ITGA2 = "#4DBBD5CC", ITGA9 = "#4DBBD5CC", ITGAV = "#4DBBD5CC", ITGB1 = "#4DBBD5CC", ITGB8 = "#4DBBD5CC", CD44 = "#4DBBD5CC", ITGA6 = "#4DBBD5CC", ITGB6 = "#4DBBD5CC",
                                         BMPR1B = "#4DBBD5CC", IGF1R = "#4DBBD5CC", PLXND1 = "#4DBBD5CC", ROBO2 = "#4DBBD5CC", 
                                         CADM1 = "#4DBBD5CC", NECTIN3 = "#4DBBD5CC", NOTCH1 = "#4DBBD5CC", NOTCH2 = "#4DBBD5CC", NRCAM = "#4DBBD5CC", CLSTN1 = "#4DBBD5CC", NLGN1 = "#4DBBD5CC")

CellChat_receptor_stroma.markers = c(ITGA1 = "#E64B35CC", ITGA2 = "#E64B35CC", ITGA9 = "#E64B35CC", ITGAV = "#E64B35CC", ITGB1 = "#E64B35CC", ITGB8 = "#E64B35CC", CD44 = "#E64B35CC", ITGA6 = "#E64B35CC", ITGB6 = "#E64B35CC",
                                     BMPR1B = "#E64B35CC", IGF1R = "#E64B35CC", PLXND1 = "#E64B35CC", ROBO2 = "#E64B35CC", 
                                     CADM1 = "#E64B35CC", NECTIN3 = "#E64B35CC", NOTCH1 = "#E64B35CC", NOTCH2 = "#E64B35CC", NRCAM = "#E64B35CC", CLSTN1 = "#E64B35CC", NLGN1 = "#E64B35CC")

# Plotting ligands in epithelium and stromal labelled bins
Epithelium_ligand_08_PCOS = Plot_region(seurat.x = seurat.list[2], cord.x = c(17000, 19000), cord.y = c(13000, 15000), marker.plot = CellChat_ligand_epithelium.markers, output.name = "CellChat_ligand_epithelium", ct.sample = "Epithelium")
Epithelium_ligand_20_LS = Plot_region(seurat.x = seurat.list[3], cord.x = c(9000, 11000), cord.y = c(11000, 13000), marker.plot = CellChat_ligand_epithelium.markers, output.name = "CellChat_ligand_epithelium", ct.sample = "Epithelium")
Epithelium_ligand_11_PCOS = Plot_region(seurat.x = seurat.list[4], cord.x = c(14000, 16000), cord.y = c(12000, 14000), marker.plot = CellChat_ligand_epithelium.markers, output.name = "CellChat_ligand_epithelium", ct.sample = "Epithelium")
Epithelium_ligand_22_Met = Plot_region(seurat.x = seurat.list[5], cord.x = c(15000, 17000), cord.y = c(14000, 16000), marker.plot = CellChat_ligand_epithelium.markers, output.name = "CellChat_ligand_epithelium", ct.sample = "Epithelium")
Epithelium_ligand_02_Ctrl = Plot_region(seurat.x = seurat.list[6], cord.x = c(11500, 13500), cord.y = c(9500, 11500), marker.plot = CellChat_ligand_epithelium.markers, output.name = "CellChat_ligand_epithelium", ct.sample = "Epithelium")
Epithelium_ligand_28_Ctrl = Plot_region(seurat.x = seurat.list[7], cord.x = c(12000, 14000), cord.y = c(8000, 10000), marker.plot = CellChat_ligand_epithelium.markers, output.name = "CellChat_ligand_epithelium", ct.sample = "Epithelium")

Stroma_ligand_08_PCOS = Plot_region(seurat.x = seurat.list[2], cord.x = c(17000, 19000), cord.y = c(13000, 15000), marker.plot = CellChat_ligand_stroma.markers, output.name = "CellChat_ligand_stromal", ct.sample = "Stromal")
Stroma_ligand_20_LS = Plot_region(seurat.x = seurat.list[3], cord.x = c(9000, 11000), cord.y = c(11000, 13000), marker.plot = CellChat_ligand_stroma.markers, output.name = "CellChat_ligand_stromal", ct.sample = "Stromal")
Stroma_ligand_11_PCOS = Plot_region(seurat.x = seurat.list[4], cord.x = c(14000, 16000), cord.y = c(12000, 14000), marker.plot = CellChat_ligand_stroma.markers, output.name = "CellChat_ligand_stromal", ct.sample = "Stromal")
Stroma_ligand_22_Met = Plot_region(seurat.x = seurat.list[5], cord.x = c(15000, 17000), cord.y = c(14000, 16000), marker.plot = CellChat_ligand_stroma.markers, output.name = "CellChat_ligand_stromal", ct.sample = "Stromal")
Stroma_ligand_02_Ctrl = Plot_region(seurat.x = seurat.list[6], cord.x = c(11500, 13500), cord.y = c(9500, 11500), marker.plot = CellChat_ligand_stroma.markers, output.name = "CellChat_ligand_stromal", ct.sample = "Stromal")
Stroma_ligand_28_Ctrl = Plot_region(seurat.x = seurat.list[7], cord.x = c(12000, 14000), cord.y = c(8000, 10000), marker.plot = CellChat_ligand_stroma.markers, output.name = "CellChat_ligand_stromal", ct.sample = "Stromal")

# Plotting receptors in epithelium and stromal labelled bins
Epithelium_receptor_08_PCOS = Plot_region(seurat.x = seurat.list[2], cord.x = c(17000, 19000), cord.y = c(13000, 15000), marker.plot = CellChat_receptor_epithelium.markers, output.name = "CellChat_receptor_epithelium", ct.sample = "Epithelium")
Epithelium_receptor_20_LS = Plot_region(seurat.x = seurat.list[3], cord.x = c(9000, 11000), cord.y = c(11000, 13000), marker.plot = CellChat_receptor_epithelium.markers, output.name = "CellChat_receptor_epithelium", ct.sample = "Epithelium")
Epithelium_receptor_11_PCOS = Plot_region(seurat.x = seurat.list[4], cord.x = c(14000, 16000), cord.y = c(12000, 14000), marker.plot = CellChat_receptor_epithelium.markers, output.name = "CellChat_receptor_epithelium", ct.sample = "Epithelium")
Epithelium_receptor_22_Met = Plot_region(seurat.x = seurat.list[5], cord.x = c(15000, 17000), cord.y = c(14000, 16000), marker.plot = CellChat_receptor_epithelium.markers, output.name = "CellChat_receptor_epithelium", ct.sample = "Epithelium")
Epithelium_receptor_02_Ctrl = Plot_region(seurat.x = seurat.list[6], cord.x = c(11500, 13500), cord.y = c(9500, 11500), marker.plot = CellChat_receptor_epithelium.markers, output.name = "CellChat_receptor_epithelium", ct.sample = "Epithelium")
Epithelium_receptor_28_Ctrl = Plot_region(seurat.x = seurat.list[7], cord.x = c(12000, 14000), cord.y = c(8000, 10000), marker.plot = CellChat_receptor_epithelium.markers, output.name = "CellChat_receptor_epithelium", ct.sample = "Epithelium")

Epithelium_receptor_08_PCOS = Plot_region(seurat.x = seurat.list[2], cord.x = c(17000, 19000), cord.y = c(13000, 15000), marker.plot = CellChat_receptor_stroma.markers, output.name = "CellChat_receptor_stromal", ct.sample = "Stromal")
Epithelium_receptor_20_LS = Plot_region(seurat.x = seurat.list[3], cord.x = c(9000, 11000), cord.y = c(11000, 13000), marker.plot = CellChat_receptor_stroma.markers, output.name = "CellChat_receptor_stromal", ct.sample = "Stromal")
Epithelium_receptor_11_PCOS = Plot_region(seurat.x = seurat.list[4], cord.x = c(14000, 16000), cord.y = c(12000, 14000), marker.plot = CellChat_receptor_stroma.markers, output.name = "CellChat_receptor_stromal", ct.sample = "Stromal")
Epithelium_receptor_22_Met = Plot_region(seurat.x = seurat.list[5], cord.x = c(15000, 17000), cord.y = c(14000, 16000), marker.plot = CellChat_receptor_stroma.markers, output.name = "CellChat_receptor_stromal", ct.sample = "Stromal")
Epithelium_receptor_02_Ctrl = Plot_region(seurat.x = seurat.list[6], cord.x = c(11500, 13500), cord.y = c(9500, 11500), marker.plot = CellChat_receptor_stroma.markers, output.name = "CellChat_receptor_stromal", ct.sample = "Stromal")
Epithelium_receptor_28_Ctrl = Plot_region(seurat.x = seurat.list[7], cord.x = c(12000, 14000), cord.y = c(8000, 10000), marker.plot = CellChat_receptor_stroma.markers, output.name = "CellChat_receptor_stromal", ct.sample = "Stromal")

# Function to overlap featureplot, plotting two markers in the same featureplot
Plot.overlay <- function(seurat.x, cord.x, cord.y, marker.plot, output.name,  
                         ct.order = celltype.order, ct.col = celltype.color, 
                         plot.cols = c("#FFFFCC", "#CC0000", "#3399FF"), 
                         seurat.ident = Seurat.label, set.pt.size = 1, set.blend = 0.25, 
                         do.subset = TRUE) {
  
  # Get the sample name
  seurat.sample = gsub(".*Output/.*//(.*?)/.*", "\\1", seurat.x)
  print(seurat.sample)
  
  # Generate output directory
  if (dir.exists(path = paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name)) == FALSE) {
    print(paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name))
    dir.create(path =paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name), recursive = TRUE)
    Output.dir = paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name, "/")
  } else if (dir.exists(path = paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name)) == TRUE) {
    print("Directory exists")
    Output.dir = paste0("Output/4_Main_plotting/", seurat.sample, "_", output.name, "/")
  } else {
    print("Error with output directory")
  }
  
  # Load the seurat object
  seurat.x = LoadH5Seurat(seurat.x)
  
  # Set the ident with selected label
  Idents(seurat.x) = seurat.ident
  
  # Reorder seurat based on selected order
  Idents(seurat.x) <- factor(Idents(seurat.x), levels= ct.order)
  
  levels(seurat.x)
  ct.col = ct.col[names(ct.col) %in% levels(seurat.x)]
  
  # Illustrade subset area
  # Do a dimplot with the clusters
  DimPlot(seurat.x, cols = ct.col, reduction = "spatial", pt.size = 1) + 
    geom_hline(yintercept=c(cord.y), linetype='dashed', col = 'black') +
    geom_vline(xintercept=c(cord.x), linetype='dashed', col = 'black') +
    theme_cowplot()
  ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_Dimplot_main_clusters_dashed_region.pdf"))
  
  # Subset the spatial dimension
  seurat.subset = subset(seurat.x, subset = spatial_1 > cord.x[1] & spatial_1 < cord.x[2] & 
                      spatial_2 > cord.y[1] & spatial_2 < cord.y[2])
  
  # Plot the marker on the full spatial and subset region
  for (marker.x in names(marker.plot)) {
    
    print(paste("Plotting", marker.x))
    
    FeaturePlot(seurat.x, features = marker.plot[[marker.x]], blend = TRUE, blend.threshold = set.blend, pt.size = 0.1, reduction = "spatial",
                cols = plot.cols)
    ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_", marker.x, "_Overlay_Featureplot_full_spatial.pdf"), dpi = 300, width = 20, height = 6)
    
    FeaturePlot(seurat.subset, features = marker.plot[[marker.x]], blend = TRUE, blend.threshold = set.blend, pt.size = set.pt.size, reduction = "spatial",
                cols = plot.cols)
    ggsave2(paste0(Output.dir,seurat.sample, "_", output.name, "_", marker.x, "_Overlay_Featureplot_selected_region.pdf"), dpi = 300, width = 20, height = 6)
    
  }
  
  
}

CellChat_ligand_stroma.markers = c(COL1A1 = "#E64B35CC", COL1A2 = "#E64B35CC", COL4A1 = "#E64B35CC", COL4A2 = "#E64B35CC", COL6A3 = "#E64B35CC", FN1 = "#E64B35CC", LAMA2 = "#E64B35CC", LAMC1 = "#E64B35CC", LAMC2 = "#E64B35CC", LAMC3 = "#E64B35CC",
                                   SPP1 = "#E64B35CC", GDF7 = "#E64B35CC", IGF1 = "#E64B35CC", SEMA3E = "#E64B35CC", SLIT2 = "#E64B35CC",
                                   CADM1 = "#E64B35CC", CNTN1 = "#E64B35CC", NRXN3 = "#E64B35CC")

# Receptor markers
CellChat_paired_markers = list(EPCAM_IGF1 = c("IGF1", "EPCAM"), 
                               COL1A2_ITGA2 = c("COL1A2", "ITGA2"),
                               CNTN1_NRCAM = c("CNTN1", "NRCAM"),
                               IGF1_IGF1R = c("IGF1", "IGF1R"))

CellChat_ECM_interactions = list(COL1A1_ITGA1 = c("COL1A1", "ITGA1"),
                                 COL1A1_ITGA2 = c("COL1A1", "ITGA2"),
                                 COL1A1_ITGB1 = c("COL1A1", "ITGB1"),
                                 COL1A1_ITGB8 = c("COL1A1", "ITGB8"),
                                 COL1A1_CD44 = c("COL1A1", "ITGB8"),
                                 COL1A2_ITGA1 = c("COL1A2", "ITGA1"),
                                 COL1A2_ITGA2 = c("COL1A2", "ITGA2"),
                                 COL1A2_ITGB1 = c("COL1A2", "ITGB1"),
                                 COL1A2_ITGB8 = c("COL1A2", "ITGB8"),
                                 COL1A2_CD44 = c("COL1A2", "ITGB8"),
                                 COL4A2_ITGA1 = c("COL4A2", "ITGA1"),
                                 COL4A2_ITGA2 = c("COL4A2", "ITGA2"),
                                 COL4A2_ITGB1 = c("COL4A2", "ITGB1"),
                                 COL4A2_ITGB8 = c("COL4A2", "ITGB8"),
                                 COL4A2_CD44 = c("COL4A2", "ITGB8"),
                                 COL6A3_ITGA1 = c("COL6A3", "ITGA1"),
                                 COL6A3_ITGA2 = c("COL6A3", "ITGA2"),
                                 COL6A3_ITGB1 = c("COL6A3", "ITGB1"),
                                 COL6A3_ITGB8 = c("COL6A3", "ITGB8"),
                                 COL6A3_CD44 = c("COL6A3", "ITGB8"),
                                 FN1_CD44 = c("FN1", "ITGAV"),
                                 FN1_CD44 = c("FN1", "ITGA6"),
                                 FN1_CD44 = c("FN1", "ITGB6"),
                                 LAMA2_ITGA6 = c("LAMA2", "ITGA6"),
                                 LAMC1_ITGA6 = c("LAMC1", "ITGA6"),
                                 LAMC2_ITGA6 = c("LAMC2", "ITGA6"),
                                 LAMC3_ITGA6 = c("LAMC3", "ITGA6"))

CellChat_Secreted_signaling = list(SPP1_ITGA9 = c("SPP1", "ITGA9"), 
                                   SPP1_ITGAV = c("SPP1", "ITGAV"),
                                   SPP1_ITGB1 = c("SPP1", "ITGB1"),
                                   IGF1_IGF1R = c("IGF1", "IGF1R"),
                                   GDF7_BMPR1B = c("GDF7", "BMPR1B"),
                                   SLIT2_ROBO2 = c("SLIT2", "ROBO2"),
                                   SEMA3E_PLXND1 = c("SEMA3E", "PLXND1"))

CellChat_Cell_Adhesion = list(CADM1_NECTIN3 = c("CADM1", "NECTIN3"), 
                                   CNTN1_NOTCH1 = c("CNTN1", "NOTCH1"),
                                   CNTN1_NOTCH2 = c("CNTN1", "NOTCH2"),
                                   CNTN1_NRCAM = c("CNTN1", "NRCAM"),
                                   NRXN3_CLSTN1 = c("NRXN3", "CLSTN1"), 
                                   NRXN3_NLGN1 = c("NRXN3", "NLGN1"))

# Plot 

# ECM repector matrix interaction
ECM_col = c("#FFFFCC", "#B2182B", "#2166AC")
sample_08_PCOS = Plot.overlay(seurat.x = seurat.list[2], cord.x = c(17000, 19000), cord.y = c(13000, 15000), marker.plot = CellChat_ECM_interactions, plot.cols = ECM_col, output.name = "ECM_interactions")
sample_027_W16 = Plot.overlay(seurat.x = seurat.list[3], cord.x = c(9000, 11000), cord.y = c(11000, 13000), marker.plot = CellChat_ECM_interactions, plot.cols = ECM_col, output.name = "ECM_interactions")
sample_11_PCOS = Plot.overlay(seurat.x = seurat.list[4], cord.x = c(14000, 16000), cord.y = c(12000, 14000), marker.plot = CellChat_ECM_interactions, plot.cols = ECM_col, output.name = "ECM_interactions")
sample_22_Met = Plot.overlay(seurat.x = seurat.list[5], cord.x = c(15000, 17000), cord.y = c(14000, 16000), marker.plot = CellChat_ECM_interactions, plot.cols = ECM_col, output.name = "ECM_interactions")
sample_02_Ctrl = Plot.overlay(seurat.x = seurat.list[6], cord.x = c(11500, 13500), cord.y = c(9500, 11500), marker.plot = CellChat_ECM_interactions, plot.cols = ECM_col, output.name = "ECM_interactions")
sample_28_Ctrl = Plot.overlay(seurat.x = seurat.list[7], cord.x = c(12000, 14000), cord.y = c(8000, 10000), marker.plot = CellChat_ECM_interactions, plot.cols = ECM_col, output.name = "ECM_interactions")

# Secreted signaling interaction
Secreted_col = c("#FFFFCC", "#C51B7D", "#4D9221")
sample_08_PCOS = Plot.overlay(seurat.x = seurat.list[2], cord.x = c(17000, 19000), cord.y = c(13000, 15000), marker.plot = CellChat_Secreted_signaling, plot.cols = Secreted_col, output.name = "Secreted_signaling")
sample_20_LS = Plot.overlay(seurat.x = seurat.list[3], cord.x = c(9000, 11000), cord.y = c(11000, 13000), marker.plot = CellChat_Secreted_signaling, plot.cols = Secreted_col, output.name = "Secreted_signaling")
sample_11_PCOS = Plot.overlay(seurat.x = seurat.list[4], cord.x = c(14000, 16000), cord.y = c(12000, 14000), marker.plot = CellChat_Secreted_signaling, plot.cols = Secreted_col, output.name = "Secreted_signaling")
sample_22_Met = Plot.overlay(seurat.x = seurat.list[5], cord.x = c(15000, 17000), cord.y = c(14000, 16000), marker.plot = CellChat_Secreted_signaling, plot.cols = Secreted_col, output.name = "Secreted_signaling")
sample_02_Ctrl = Plot.overlay(seurat.x = seurat.list[6], cord.x = c(11500, 13500), cord.y = c(9500, 11500), marker.plot = CellChat_Secreted_signaling, plot.cols = Secreted_col, output.name = "Secreted_signaling")
sample_28_Ctrl = Plot.overlay(seurat.x = seurat.list[7], cord.x = c(12000, 14000), cord.y = c(8000, 10000), marker.plot = CellChat_Secreted_signaling, plot.cols = Secreted_col, output.name = "Secreted_signaling")

# Cell adhesion molecules interactions
Adhesion_col = c("#FFFFCC", "#B35806", "#542788")
sample_08_PCOS = Plot.overlay(seurat.x = seurat.list[2], cord.x = c(17000, 19000), cord.y = c(13000, 15000), marker.plot = CellChat_Cell_Adhesion, plot.cols = Adhesion_col, output.name = "Cell_Adhesion")
sample_20_LS = Plot.overlay(seurat.x = seurat.list[3], cord.x = c(9000, 11000), cord.y = c(11000, 13000), marker.plot = CellChat_Cell_Adhesion, plot.cols = Adhesion_col, output.name = "Cell_Adhesion")
sample_11_PCOS = Plot.overlay(seurat.x = seurat.list[4], cord.x = c(14000, 16000), cord.y = c(12000, 14000), marker.plot = CellChat_Cell_Adhesion, plot.cols = Adhesion_col, output.name = "Cell_Adhesion")
sample_22_Met = Plot.overlay(seurat.x = seurat.list[5], cord.x = c(15000, 17000), cord.y = c(14000, 16000), marker.plot = CellChat_Cell_Adhesion, plot.cols = Adhesion_col, output.name = "Cell_Adhesion")
sample_02_Ctrl = Plot.overlay(seurat.x = seurat.list[6], cord.x = c(11500, 13500), cord.y = c(9500, 11500), marker.plot = CellChat_Cell_Adhesion, plot.cols = Adhesion_col, output.name = "Cell_Adhesion")
sample_28_Ctrl = Plot.overlay(seurat.x = seurat.list[7], cord.x = c(12000, 14000), cord.y = c(8000, 10000), marker.plot = CellChat_Cell_Adhesion, plot.cols = Adhesion_col, output.name = "Cell_Adhesion")

