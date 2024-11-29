#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
#library(celldex)
library(openxlsx)


# Setting up script
Project_name = "Endo_All_Epithelium"
Input.dir = "Output/5_Epithelium_Clustering/"
select.idents = "Epithelium_labelled"

#Project_name = "Endo_All_Stromal_uSMC"
#Input.dir = "Output/5_Stromal_uSMC_Clustering/"
#select.idents = "Stromal_labelled"

#Project_name = "Endo_All_Immune"
#Input.dir = "Output/5_Immune_Clustering/"
#select.idents = "Immune_labelled"

#Project_name = "Endo_All_Endothelial"
#Input.dir = "Output/5_Endothelial_Clustering/"
#select.idents = "Endothelial_labelled"

theme_set(theme_cowplot())

FindMarker.flag = TRUE
save.subset.seurat = FALSE
MAST.test.flag = TRUE # Default is TRUE

# Checking and/or generating output dir
if (dir.exists(path = paste0("Output/6_DGE_analysis_", select.idents)) == FALSE) {
  print(paste0("Generating output directory Output/6_DGE_analysis_", select.idents))
  dir.create(path = paste0("Output/6_DGE_analysis_", select.idents), recursive = TRUE)
  Output.dir = paste0("Output/6_DGE_analysis_", select.idents)
} else if (dir.exists(path = paste0("Output/6_DGE_analysis_", select.idents)) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/6_DGE_analysis_", select.idents)
} else {
  print("Error with output directory")
}

# Loading reclustered Seurat
print("Seurat object loading")
endo.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_reclustered_labelled.h5seurat"))

# Re-order the object for plotting
# Hash out for testing if it is correct
#Idents(subset.cell) <- "Group_Stage"
#subset.cell$Group_Stage <- factor(subset.cell$Group_Stage, levels= treatment.levels)
patient.levels = c("1_Control", "2_Control", "3_Control", "4_Control", "5_Control",
                   "6_PCOS_Lifestyle_W0", "7_PCOS_Lifestyle_W0", "8_PCOS_Lifestyle_W0", "9_PCOS_EA_W0",
                   "10_PCOS_Metformin_W0", "11_PCOS_Metformin_W0",  "11_PCOS_Metformin_W0",
                   "12_PCOS_Metformin_W0",  "13_PCOS_Metformin_W0", "14_PCOS_Metformin_W0",
                   "15_PCOS_Metformin_W0", "16_PCOS_Metformin_W0",  "17_PCOS_Metformin_W0",
                   "18_PCOS_Lifestyle_W16", "19_PCOS_Lifestyle_W16", "20_PCOS_Lifestyle_W16",
                   "21_PCOS_Metformin_W16", "22_PCOS_Metformin_W16", "23_PCOS_Metformin_W16",
                   "24_PCOS_Metformin_W16", "25_PCOS_Metformin_W16", "26_PCOS_Metformin_W16",
                   "27_PCOS_Metformin_W16")
endo.integrated$orig.ident <- factor(endo.integrated$orig.ident, levels= patient.levels)

# Subset only the Control and PCOS W0. Redo normalisation, scaling, pca and UMAP.
Idents(object = endo.integrated) <- "Group_Stage"
DefaultAssay(endo.integrated) = "RNA"

endo.integrated.Ctrl_PCOS = subset(endo.integrated, idents = c("Control", "PCOS_W0"))
endo.integrated.PCOS_Metformin = subset(endo.integrated, idents = c("PCOS_W0", "PCOS_W16_Met"))
endo.integrated.PCOS_Lifestyle = subset(endo.integrated, idents = c("PCOS_W0", "PCOS_W16_LS"))

# Setting default setting for the seurat objects
#Idents(object = endo.integrated) <- select.idents
Idents(object = endo.integrated.Ctrl_PCOS) <- select.idents
Idents(object = endo.integrated.PCOS_Metformin) <- select.idents
Idents(object = endo.integrated.PCOS_Lifestyle) <- select.idents
Norm.assay = "RNA"

treatment.levels.CtrlvsPCOS = c("Control", "PCOS_W0")
x.labels.CtrlvsPCOS = c("Control", "PCOS")

treatment.levels.PCOSvsMetformin = c("PCOS_W16_Met", "PCOS_W0")
x.labels.PCOSvsMetformin = c("Metformin", "PCOS")

treatment.levels.PCOSvsLifestyle = c("PCOS_W16_LS", "PCOS_W0")
x.labels.PCOSvsLifestyle = c("Lifestyle", "PCOS")

###### EDITING FUNCTION #########
Run_FindMarkers_Specific <- function(seurat.x = NULL, celltype, x = "PCOS_W0", y = "Control",
                                   group = "Group_Stage", cluster = select.idents, min.cell.exp = 0.25,
                                   MAST.test = TRUE, selected.genes = NULL, Project_name = "Endo_All_CtrlvsPCOS",
                                   treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                   x.labels = c("Control", "PCOS", "Metformin", "Lifestyle")) {
  

  if (dir.exists(path = paste0(Output.dir, "/", celltype)) == FALSE) {
    print(paste0("Generating output directory"))
    dir.create(path = paste0(Output.dir, "/", celltype), recursive = TRUE)
    Output.dir.celltype = paste0(Output.dir, "/", celltype, "/")
  } else if (dir.exists(path = paste0(Output.dir, "/", celltype)) == TRUE) {
    print("Directory exists")
    Output.dir.celltype = paste0(Output.dir, "/", celltype, "/")
  } else {
    print("Error with output directory")
  }
  
  # Setting plotting colors
  Heatmap.cols = c("dodgerblue", "seashell", "firebrick")
  Vln.cols = c("dodgerblue", "firebrick", "forestgreen", "sandybrown")
  Vln.idents.cols = c(rep("dodgerblue", 5), rep("firebrick", 12), 
                      rep("forestgreen", 3), rep("sandybrown", 7))
  Dot.cols = c("dodgerblue", "firebrick")
  
  # Subsetting target celltype
  subset.cell = subset(seurat.x, idents = celltype)
  
  # NULL number of DEGs
  # Print total of DEG's
  n.DEG = NULL
  n.DEG.sign = NULL
  
  #subset.cell.subsampled <- subset.cell[, sample(colnames(subset.cell), size =1000, replace=F)]
  
  # Rescale data after subsetting
  # Normalise and rescale data after subsetting
  print("Setting DefaultAssay to RNA and log2 normalising it")
  DefaultAssay(subset.cell) = "RNA"
  Idents(subset.cell) <- "orig.ident"
  #Idents(object = subset.cell) <- "orig.ident" # TESTING NORMALIZING OVER IDENTS
  subset.cell = NormalizeData(subset.cell)
  subset.cell = FindVariableFeatures(subset.cell)
  all.genes = rownames(subset.cell)
  subset.cell = ScaleData(subset.cell, features = all.genes)  #endo.integrated = FindVariableFeatures(endo.integrated)
  
  # Run FindMarkers
  Idents(subset.cell) <- "Group_Stage"
  
  if (MAST.test == TRUE) {
    print("FindMarker with MAST hurdle model")
    
    DEG_list_x.y = FindMarkers(object = subset.cell, ident.1 = x, ident.2 = y,
                               group.by = "Group_Stage", min.pct = min.cell.exp,
                               test.use = "MAST")
    
  } else if (MAST.test == FALSE) {
    print("FindMarker with Wilcoxon Rank Sum test")
    
    DEG_list_x.y = FindMarkers(object = subset.cell, ident.1 = x, ident.2 = y,
                               group.by = "Group_Stage", min.pct = min.cell.exp,
                               test.use = "wilcox")
  }
  
  n.DEG = nrow(DEG_list_x.y)
  
  # Extract significant genes
  DEG_sign_x.y = DEG_list_x.y[DEG_list_x.y$p_val_adj < 0.05,]
  
  if (nrow(DEG_sign_x.y) == 0) {
    print("No significant DEGs")
    write.xlsx(DEG_list_x.y, 
               paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_", n.DEG, "DEGs_all_DEG.xlsx"),
               colNames = TRUE, rowNames = TRUE)
    return(NA)
  }
  
  n.DEG.sign = nrow(DEG_sign_x.y)
  
  # Save the DEGs as a rds
  write.xlsx(DEG_list_x.y, 
             paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_", n.DEG, "DEGs_all_DEG.xlsx"),
             colNames = TRUE, rowNames = TRUE)
  write.xlsx(DEG_sign_x.y, 
             paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_", n.DEG.sign, "DEGs_sign_DEG.xlsx"),
             colNames = TRUE, rowNames = TRUE)
  
  # Order DEG list based on log2FC
  DEG_sign_x.y = DEG_sign_x.y[order(DEG_sign_x.y$avg_log2FC),]

  # Extract DEG names
  DEG_ID_x.y = rownames(DEG_sign_x.y)

  # Plotting when only control and one treatment is compared to baseline
  
  if (nrow(DEG_sign_x.y) >= 10) {
    print("More than 10 DEG's")
    
    # Generating a dotplot of all DEG's
    plot.res = DotPlot(subset.cell, features = rownames(DEG_sign_x.y), cols = Dot.cols ,dot.scale = 12, scale = FALSE) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) + theme(axis.text.y = element_text(face = "italic")) +
      RotatedAxis() + coord_flip() # + scale_y_discrete(labels = x.labels)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Dotplot_Specific_All.pdf"), dpi = 700)
    
    # Generating a heatmap of all DEG's
    plot.All.DEG.up = DEG_sign_x.y[DEG_sign_x.y$avg_log2FC > 0,]
    plot.All.DEG.up = plot.All.DEG.up[order(plot.All.DEG.up$avg_log2FC, decreasing = TRUE),]
    plot.All.DEG.down = DEG_sign_x.y[DEG_sign_x.y$avg_log2FC < 0,]
    plot.All.DEG.down = plot.All.DEG.down[order(plot.All.DEG.down$avg_log2FC, decreasing = TRUE),]
    
    plot.res = DoHeatmap(subset.cell, features = c(rownames(plot.All.DEG.up), rownames(plot.All.DEG.down)), 
                         raster = FALSE, group.bar = TRUE, draw.lines = TRUE, 
                         angle = 0, hjust = 0.5, group.bar.height = 0.01) + 
      scale_fill_gradientn(colors = c("dodgerblue", "seashell", "firebrick")) + theme(axis.text.y = element_text(face = "italic"))
    
    ggsave2(plot = plot.res, paste0(Output.dir.celltype, Project_name,"_", Norm.assay, "_", celltype, "_All_DEG_Heatmap.pdf"), 
            dpi = 700, height = 14, width = 8)
    
    plot.res = DoHeatmap(subset.cell, features = c(rownames(plot.All.DEG.up), rownames(plot.All.DEG.down)), 
                         group.by = "orig.ident", raster = FALSE, group.bar = TRUE, draw.lines = TRUE, 
                         angle = 0, hjust = 0.5, group.bar.height = 0.01) + 
      scale_fill_gradientn(colors = c("dodgerblue", "seashell", "firebrick")) + theme(axis.text.y = element_text(face = "italic"))
    
    ggsave2(plot = plot.res, paste0(Output.dir.celltype, Project_name,"_", Norm.assay, "_", celltype, "_All_DEG_Heatmap_orig_idents.pdf"), 
            dpi = 700, height = 14, width = 8)
    
    # Extract top 5 and bottom 5 
    if (nrow(plot.All.DEG.up) > 5) {
      plot.top.DEG.up = plot.All.DEG.up[1:5,]
    } else if (nrow(plot.All.DEG.up) <= 5) {
      plot.top.DEG.up = plot.All.DEG.up
      plot.top.DEG.up.length = nrow(plot.top.DEG.up)
    }
    
    if (nrow(plot.All.DEG.down) > 5) {
      plot.top.DEG.down = plot.All.DEG.down[(nrow(plot.All.DEG.down)-4):nrow(plot.All.DEG.down),]
    } else if (nrow(plot.All.DEG.down) <= 5) {
      plot.top.DEG.down = plot.All.DEG.down
      plot.top.DEG.down.length = nrow(plot.top.DEG.down)
    }
    
    # Generating dotplots
    plot.res = DotPlot(subset.cell, features = c(rownames(plot.top.DEG.down), rownames(plot.top.DEG.up), scale = FALSE), 
                       cols = Dot.cols, dot.scale = 12, assay = "RNA") + theme(axis.text.y = element_text(face = "italic")) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) + RotatedAxis() + coord_flip() # +
      #scale_y_discrete(labels = x.labels)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Dotplot_Specific_top-down.pdf"), dpi = 700)
    
    plot.res = DotPlot(subset.cell, features = c(rownames(plot.top.DEG.down), rownames(plot.top.DEG.up)), cols = Dot.cols, dot.scale = 10, scale = FALSE) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_text(face = "italic")) #+
      #scale_y_discrete(labels = x.labels)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Dotplot_flipped_Specific_top-down.pdf"), dpi = 700)
    
    # Generating violin plots
    plot.res = VlnPlot(subset.cell, features = c(rownames(plot.top.DEG.down), rownames(plot.top.DEG.up)), pt.size = 0, cols = Vln.cols)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Vlnplot_Specific_top-down.pdf"), dpi = 700, 
            height = 10, width = 10)
    
    plot.res = VlnPlot(subset.cell, features = c(rownames(plot.top.DEG.down), rownames(plot.top.DEG.up)), pt.size = 0.1, cols = Vln.cols)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Vlnplot_Specific_dots_top-down.pdf"), dpi = 700, 
            height = 10, width = 10)
    
    plot.res = VlnPlot(subset.cell, features = c(rownames(plot.top.DEG.down), rownames(plot.top.DEG.up)), pt.size = 0.1, group.by = "orig.ident", cols = Vln.idents.cols)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Vlnplot_Specific_top-down_orig_ident.pdf"), dpi = 700, 
            height = 10, width = 17)
    
    # Generate heatmaps
    plot.res = DoHeatmap(subset.cell, features = c(rownames(plot.top.DEG.down), rownames(plot.top.DEG.up)), raster = FALSE,
                         group.bar = TRUE, group.bar.height = 0.01, draw.lines = TRUE, hjust = 0.5, angle = 0) + 
      scale_fill_gradientn(colors = Heatmap.cols) + theme(axis.text.y = element_text(face = "italic"))
    
    ggsave2(plot = plot.res, paste0(Output.dir.celltype, Project_name,"_", Norm.assay, "_", celltype, "_top_DEG_Heatmap.pdf"), 
            dpi = 700, height = 7, width = 6)
    
    plot.res = DoHeatmap(subset.cell, features = c(rownames(plot.top.DEG.down), rownames(plot.top.DEG.up)), raster = FALSE,
                         group.by = "orig.ident", group.bar = TRUE, group.bar.height = 0.01, draw.lines = TRUE, hjust = 0.5, angle = 0) + 
      scale_fill_gradientn(colors = Heatmap.cols) + theme(axis.text.y = element_text(face = "italic"))
    
    ggsave2(plot = plot.res, paste0(Output.dir.celltype, Project_name,"_", Norm.assay, "_", celltype, "_top_DEG_Heatmap_orig_idents.pdf"), 
            dpi = 700, height = 7, width = 6)
    
    plot.res = NULL
  } else if (nrow(DEG_sign_x.y) < 10) {
    print("Less than 10 DEG's")
    plot.DEG = rownames(DEG_sign_x.y)
    
    DotPlot(subset.cell, features = plot.DEG)
    
    # Dotplots are not scaled as low number of groups may produce misleading results. 
    # Therefore, violinplots are produced are a better way of illustrating the DEGs. 
    
    plot.res = DotPlot(subset.cell, features = plot.DEG, cols = Dot.cols, dot.scale = 12, scale = FALSE) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) + theme(axis.text.y = element_text(face = "italic")) +
      RotatedAxis() + coord_flip() + scale_y_discrete(labels = x.labels)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Dotplot_Specific_All.pdf"), dpi = 700)
    
    plot.res = DotPlot(subset.cell, features = plot.DEG, cols = Dot.cols, dot.scale = 10, scale = FALSE) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_text(face = "italic")) + scale_y_discrete(labels = x.labels)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Dotplot_flipped_Specific_All.pdf"), dpi = 700)
    
    plot.res = VlnPlot(subset.cell, features = plot.DEG, pt.size = 0, cols = Vln.cols)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Vlnplot_Specific_All.pdf"), dpi = 700, 
            height = 10, width = 10)
    
    plot.res = VlnPlot(subset.cell, features = plot.DEG, pt.size = 0.1, cols = Vln.cols)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Vlnplot_Specific_dots_All.pdf"), dpi = 700, 
            height = 10, width = 10)
    
    plot.res = VlnPlot(subset.cell, features = plot.DEG, pt.size = 0.1, group.by = "orig.ident", cols = Vln.idents.cols)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Vlnplot_Specific_All_orig_ident.pdf"), dpi = 700, 
            height = 10, width = 10)
    
    plot.res = NULL
  }
  
  # Plot selected genes
  if (!is.null(selected.genes) == TRUE) {
    
    selected.genes = unique(selected.genes)
      
    plot.res = DotPlot(subset.cell, features = selected.genes, cols = Dot.cols, dot.scale = 12, scale = FALSE) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) + theme(axis.text.y = element_text(face = "italic")) +
      RotatedAxis() + coord_flip() + scale_y_discrete(labels = x.labels)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Dotplot_Specific_All_selected_genes.pdf"), dpi = 700)
    
    plot.res = DotPlot(subset.cell, features = selected.genes, cols = Dot.cols, dot.scale = 12, scale = FALSE) +
      theme(axis.text.x=element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_text(face = "italic")) +
      scale_y_discrete(labels = x.labels)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Dotplot_flipped_Specific_All_selected_genes.pdf"), dpi = 700)
    
    plot.res = VlnPlot(subset.cell, features = selected.genes, pt.size = 0, cols = Vln.cols)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Vlnplot_Specific_All_selected_genes.pdf"), dpi = 700, 
            height = 10, width = 10)
    
    plot.res = VlnPlot(subset.cell, features = selected.genes, pt.size = 0.1, cols = Vln.cols)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Vlnplot_Specific_dots_All_selected_genes.pdf"), dpi = 700, 
            height = 10, width = 10)
    
    plot.res = VlnPlot(subset.cell, features = selected.genes, pt.size = 0.1, group.by = "orig.ident", cols = Vln.idents.cols)
    ggsave2(plot = plot.res, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_Vlnplot_Specific_All_orig_ident_selected_genes.pdf"), dpi = 700, 
            height = 10, width = 10)
    
    plot.res = NULL
  }
  
  # Print total of DEG's
  print(paste(nrow(DEG_list_x.y), "DEG's in", celltype))
  print(paste(nrow(DEG_sign_x.y), "significant DEG's in", celltype))
  
  # Save the seurat object for future analysis
  if (save.subset.seurat == TRUE) {
    SaveH5Seurat(subset.cell, paste0(Output.dir.celltype,Project_name, "_", celltype, "_", x, "_", y, "_seurat.h5seurat"),
                 overwrite = TRUE)
  }
  
}

# Run functions
if (select.idents == "Immune_labelled") {
  print("Running PCOS vs. Control")
  ### Run FindMarkers on Ctrl vs. PCOS W0 ####
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "T-cells CD4+", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x = "PCOS_W0", y = "Control",
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = c("PGR", "ESR1", "AR"))
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "T-cells CD8+", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x = "PCOS_W0", y = "Control",
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "uNK 1", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x = "PCOS_W0", y = "Control",
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "uNK 2", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "uNK 3", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels =  x.labels.CtrlvsPCOS,
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "uM 1", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "uM 2", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = c("PGR", "ESR1", "AR"))
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "B-cells", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  # No DEGs found in celltypes below
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Tregs", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "ILC3", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Mast cells", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "DC1", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "DC2", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Migratory DC", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "pDC", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Undefined", 
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           x = "PCOS_W0", y = "Control",
                           selected.genes = NULL)
  ##########
  
  ### Run FindMarkers on PCOS vs. Metformin ####
  print("Running Metformin vs. PCOS")
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "T-cells CD4+", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = c("PGR", "ESR1", "AR"))
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "T-cells CD8+", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "uNK 1", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "uNK 2", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "uNK 3", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels =  x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "uM 1", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "uM 2", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = c("PGR", "ESR1", "AR"))
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "B-cells", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Tregs", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "ILC3", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Mast cells", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "DC1", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "DC2", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Migratory DC", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "pDC", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Undefined", 
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsMetformin",
                           x.labels = x.labels.PCOSvsMetformin,
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           selected.genes = NULL)
  ##########
  
  ### Run FindMarkers on PCOS vs. Lifestyle ####
  print("Running Lifestyle vs. PCOS")
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "T-cells CD4+", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = c("PGR", "ESR1", "AR"))
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "T-cells CD8+", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "uNK 1", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "uNK 2", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "uNK 3", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           x.labels =  x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "uM 1", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "uM 2", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = c("PGR", "ESR1", "AR"))
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "B-cells", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Tregs", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "ILC3", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Mast cells", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "DC1", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "DC2", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Migratory DC", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "pDC", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Undefined", 
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           x.labels = x.labels.PCOSvsLifestyle,
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           selected.genes = NULL)
  ##########
  
} else if (select.idents == "Epithelium_labelled") {
  
  ### Run FindMarkers on Ctrl vs. PCOS W0 ####
  print("FindMarkers on Ctrl vs. PCOS")
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Lumenal",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = c("SEMA3E", "ROBO2", "NRCAM", "GPC6", "CTNNA2", "HMCN1", #Downregulated
                                              "SLPI", "ADAMTS9", "SPON1", "ITIH5", "LAMC2", "THSD4")) # Upregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "SOX9+ LGR5+",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = c("COL3A1", "ESR1", "CTNNA2", "CENPK", "CENPP", "ACSM1", "ACSM3", "BRCA1", #Downregulated
                                              "CD44", "CD55", "ITGA2", "ITGA3", "ITGB3", "PAEP", "MMP3", "MMP10")) # Upregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "SOX9+ LGR5-",
                           x = "PCOS_W0", y = "Control",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = c("SEMA3E", "CTNNA2", "AUTS2", "NRCAM", "FGF13", "GPC6", "SOX5", "ZFP36L2",
                                              "CNTN5", "UNC5C", "SFRP4", "VIM", "CDH6")) #Downregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Glandular",
                           x = "PCOS_W0", y = "Control",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels =  x.labels.CtrlvsPCOS,
                           selected.genes = c("ESR1", "LEF1", "NRCAM", "BMPR1B", "GREM2", "SOX5", "COL3A1", "GPC6", #Downregulated
                                              "SOD2", "GLI3", "HNF1B", "AIMP1", "ACSL4", "FOXO1", "SLC1A1", "SLC15A1",
                                              "TXNIP", "ITGA2", "PAX8", "SLC15A4", "ERN1", "LEPR", "MGST1", "SLC7A2", "SLC44A1", "GPX3", "MAP3K5")) #Upregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "SOX9+ proliferative",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = NULL)

  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "AR+",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = c("MMP11", "LAMA2", "ROBO2", "GPC6", "ECM2", "COL10A1", "SEMA3E", "CNTN1", "ESR1", #Downregulated
                                              "MMP3", "MMP10", "ITGA2", "SLPI", "COL7A1", "CD44", "STC1", "NEAT1")) # Upregulated

  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Ciliated",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = c("CFAP43", "CFAP58", "CFAP73", "CFAP157", "RP1", "FLACC1", #Downregulated
                                              "SVIL", "RCC1", "CDK14", "ITGA6", "MYH9", "E2F7")) # Upregulated
  ##########
  
  ### Run FindMarkers on PCOS vs. Metformin ####
  print("FindMarkers on PCOS vs. Metformin")

  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Lumenal",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = c("CTNNA2", "COL3A1", "TMLHE", "MAP2K6", "STXBP6", 
                                              "DLGAP1", "CLDN10", "NEAT1", "CD55", "PAEP"))
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "SOX9+ LGR5+",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "SOX9+ LGR5-",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Glandular",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels =  x.labels.PCOSvsMetformin,
                           selected.genes = c("HIBCH", "MFSD4B", "STX18", "ANKUB1", "ESR1", 
                                              "CD55", "MUC16", "CD44", "NEAT1", "MALAT1"))
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "SOX9+ proliferative",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "AR+",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Ciliated",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  ##########
  
  ### Run FindMarkers on PCOS vs. Lifestyle ####
  print("FindMarkers on PCOS vs. Lifestyle")
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Lumenal",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = c("CTNNA2", "COL3A1", "TMLHE", "MAP2K6", "STXBP6", 
                                              "DLGAP1", "CLDN10", "NEAT1", "CD55", "PAEP"))
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "SOX9+ LGR5+",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "SOX9+ LGR5-",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Glandular",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels =  x.labels.PCOSvsLifestyle,
                           selected.genes = c("HIBCH", "MFSD4B", "STX18", "ANKUB1", "ESR1", 
                                              "CD55", "MUC16", "CD44", "NEAT1", "MALAT1"))
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "SOX9+ proliferative",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "AR+",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Ciliated",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
} else if (select.idents == "Stromal_labelled") {
  
  ### Run FindMarkers on Ctrl vs. PCOS W0 ####
  print("FindMarkers on Ctrl vs. PCOS")
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Stroma 1",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = NULL) # Upregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Stroma 2",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = NULL) # Upregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Stroma proliferative",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = NULL) # Upregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Fibroblast",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = NULL) # Upregulated
  
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "uSMC",
                           x = "PCOS_W0", y = "Control",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = NULL) #Downregulated
  
  ##########
  
  ### Run FindMarkers on PCOS vs. Metformin ####
  print("FindMarkers on PCOS vs. Metformin")
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Stroma 1",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Stroma 2",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Stroma proliferative",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Fibroblast",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "uSMC",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
 
  ##########
  
  ### Run FindMarkers on PCOS vs. Lifestyle ####
  print("FindMarkers on PCOS vs. Lifestyle")
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Stroma 1",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Stroma 2",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Stroma proliferative",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Fibroblast",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "uSMC",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
} else if (select.idents == "Endothelial_labelled") {
  
  ### Run FindMarkers on Ctrl vs. PCOS W0 ####
  print("FindMarkers on Ctrl vs. PCOS")
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Endothelial Vein",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = c("SEMA3E", "ROBO2", "NRCAM", "GPC6", "CTNNA2", "HMCN1", #Downregulated
                                              "SLPI", "ADAMTS9", "SPON1", "ITIH5", "LAMC2", "THSD4")) # Upregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Endothelial Artery",
                           x = "PCOS_W0", y = "Control",
                           Project_name = "Endo_All_CtrlvsPCOS",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = c("COL3A1", "ESR1", "CTNNA2", "CENPK", "CENPP", "ACSM1", "ACSM3", "BRCA1", #Downregulated
                                              "CD44", "CD55", "ITGA2", "ITGA3", "ITGB3", "PAEP", "MMP3", "MMP10")) # Upregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Endothelial proliferative",
                           x = "PCOS_W0", y = "Control",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels = x.labels.CtrlvsPCOS,
                           selected.genes = c("SEMA3E", "CTNNA2", "AUTS2", "NRCAM", "FGF13", "GPC6", "SOX5", "ZFP36L2",
                                              "CNTN5", "UNC5C", "SFRP4", "VIM", "CDH6")) #Downregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Mesenchymal",
                           x = "PCOS_W0", y = "Control",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels =  x.labels.CtrlvsPCOS,
                           selected.genes = c("ESR1", "LEF1", "NRCAM", "BMPR1B", "GREM2", "SOX5", "COL3A1", "GPC6", #Downregulated
                                              "SOD2", "GLI3", "HNF1B", "AIMP1", "ACSL4", "FOXO1", "SLC1A1", "SLC15A1",
                                              "TXNIP", "ITGA2", "PAX8", "SLC15A4", "ERN1", "LEPR", "MGST1", "SLC7A2", "SLC44A1", "GPX3", "MAP3K5")) #Upregulated
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.Ctrl_PCOS, celltype = "Lymphatic",
                           x = "PCOS_W0", y = "Control",
                           treatment.levels = treatment.levels.CtrlvsPCOS, MAST.test = MAST.test.flag,
                           Project_name = "Endo_All_CtrlvsPCOS",
                           x.labels =  x.labels.CtrlvsPCOS,
                           selected.genes = c("ESR1", "LEF1", "NRCAM", "BMPR1B", "GREM2", "SOX5", "COL3A1", "GPC6", #Downregulated
                                              "SOD2", "GLI3", "HNF1B", "AIMP1", "ACSL4", "FOXO1", "SLC1A1", "SLC15A1",
                                              "TXNIP", "ITGA2", "PAX8", "SLC15A4", "ERN1", "LEPR", "MGST1", "SLC7A2", "SLC44A1", "GPX3", "MAP3K5")) #Upregulated
  
  
  ##########
  
  ### Run FindMarkers on PCOS vs. Metformin ####
  print("FindMarkers on PCOS vs. Metformin")
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Endothelial Vein",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = c("CTNNA2", "COL3A1", "TMLHE", "MAP2K6", "STXBP6", 
                                              "DLGAP1", "CLDN10", "NEAT1", "CD55", "PAEP"))
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Endothelial Artery",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Endothelial proliferative",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsMetformin,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Mesenchymal",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels =  x.labels.PCOSvsMetformin,
                           selected.genes = c("HIBCH", "MFSD4B", "STX18", "ANKUB1", "ESR1", 
                                              "CD55", "MUC16", "CD44", "NEAT1", "MALAT1"))
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Metformin, celltype = "Lymphatic",
                           x = "PCOS_W16_Met", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsMetformin",
                           treatment.levels = treatment.levels.PCOSvsMetformin, MAST.test = MAST.test.flag,
                           x.labels =  x.labels.PCOSvsMetformin,
                           selected.genes = c("HIBCH", "MFSD4B", "STX18", "ANKUB1", "ESR1", 
                                              "CD55", "MUC16", "CD44", "NEAT1", "MALAT1"))
  
  ##########
  
  ### Run FindMarkers on PCOS vs. Lifestyle ####
  print("FindMarkers on PCOS vs. Lifestyle")
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Endothelial Vein",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = c("CTNNA2", "COL3A1", "TMLHE", "MAP2K6", "STXBP6", 
                                              "DLGAP1", "CLDN10", "NEAT1", "CD55", "PAEP"))
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Endothelial Artery",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Endothelial proliferative",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels = x.labels.PCOSvsLifestyle,
                           selected.genes = NULL)
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Mesenchymal",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels =  x.labels.PCOSvsLifestyle,
                           selected.genes = c("HIBCH", "MFSD4B", "STX18", "ANKUB1", "ESR1", 
                                              "CD55", "MUC16", "CD44", "NEAT1", "MALAT1"))
  
  Run_FindMarkers_Specific(seurat.x = endo.integrated.PCOS_Lifestyle, celltype = "Lymphatic",
                           x = "PCOS_W16_LS", y = "PCOS_W0",
                           Project_name = "Endo_All_PCOSvsLifestyle",
                           treatment.levels = treatment.levels.PCOSvsLifestyle, MAST.test = MAST.test.flag,
                           x.labels =  x.labels.PCOSvsLifestyle,
                           selected.genes = c("HIBCH", "MFSD4B", "STX18", "ANKUB1", "ESR1", 
                                              "CD55", "MUC16", "CD44", "NEAT1", "MALAT1"))
  
}



if (FindMarker.flag == TRUE) {
  print("Finding markers")
  Idents(object = endo.integrated) <- "Epithelium_labelled" #Set to analyse clusters
  endo.markers <-FindAllMarkers(endo.integrated, assay = "RNA")
  print("Saving markers as .csv")
  write.csv(endo.markers, paste0(Output.dir, Project_name, "_FindAllMarkers.csv"), quote = F)
  print("Saving markers as .rds")
  saveRDS(endo.markers, paste0(Output.dir, Project_name, "_Endo_markers.rds"))
  print("Extracting top 10 markers per cluster")
  top10 = endo.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
  
  
  print("Generating heatmap")
  endo.heatmap = DoHeatmap(endo.integrated, features = top10$gene) + NoLegend()
  ggsave2(filename = paste0(Output.dir, Project_name, "_Top10_Genes_Heatmap.pdf"),
          plot = endo.heatmap,
          dpi = 700)
} 
 
if (FindMarker.flag == TRUE) {
  print("Finding markers")
  Idents(object = endo.integrated) <- select.idents #Set to analyse clusters
  endo.markers <-FindAllMarkers(endo.integrated, assay = "RNA")
  print("Saving markers as .csv")
  write.csv(endo.markers, paste0(Output.dir, Project_name, "_FindAllMarkers.csv"), quote = F)
  print("Saving markers as .rds")
  saveRDS(endo.markers, paste0(Output.dir, Project_name, "_Endo_markers.rds"))
  print("Extracting top 10 markers per cluster")
  top10 = endo.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
  
  
  print("Generating heatmap")
  endo.heatmap = DoHeatmap(endo.integrated, features = top10$gene) + NoLegend()
  ggsave2(filename = paste0(Output.dir, Project_name, "_Top10_Genes_Heatmap.pdf"),
          plot = endo.heatmap,
          dpi = 700)
}
print("Done")
