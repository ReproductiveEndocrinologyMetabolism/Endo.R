#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(stringr)
library(openxlsx)
library(pheatmap)

#celltype = "Stromal"
plot.celltype = "Immune" # Epithelium_All, Immune, Endothelial, Epithelium_CtrlvsPCOS, Stromal_uSMC
#plot.subtype = c("AR+", "SOX9+ LGR5+", "SOX9+ proliferative", "Fibroblast", "Stroma proliferative", "uSMC")
#lot.subtype = c("Endothelial Artery", "Lymphatic")
plot.subtype = c("uNK 1", "uNK 2", "uNK 3", "uM 1", "uM 2")
#plot.genes = c("ESR1", "PGR", "AR")
#plot.genes = "ESR1"
#plot.genes = "PGR"
#plot.genes = c("PDE3A", "MGP", "FN1", "SULF1", "CLDN5", "FTL", "HDAC9", "RCAN2")
#plot.genes = c("PDE3A", "CLDN5", "FN1", "SULF1", "CD9", "DSCAM")
plot.genes = "PAEP"
Project_name = paste0("Endo_Baseline_Immune_PAEP_")
seurat.dir = "Data/Seurat_Subset/"
seurat.ident = "Group_Stage"
Group.cols = c("#A0A0A0", "#D098B2", "#95BFE1", "#65D46E")
do.boxplot = FALSE
do.Ctrl_PCOS = FALSE
#Project_name = paste0("Endo_All_", celltype, "_All")
#GO.dir = "Data/Selected_GO/"
#DEG.dir = "Data/DEG_tables/"
Norm.assay = "RNA"
Group.order = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS")

if (dir.exists(path = paste0("Output/8_Violinplot")) == FALSE) {
  print(paste0("Generating output directory Output/8_Violinplot"))
  dir.create(path = paste0("Output/8_Violinplot"), recursive = TRUE)
  Output.dir = paste0("Output/8_Violinplot")
} else if (dir.exists(path = paste0("Output/8_Violinplot")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/8_Violinplot/")
} else {
  print("Error with output directory")
}

# Loop for generating violinplots of selected genes in selected subtypes
Vln.all.list = list()
Vln.all.gene.list = list()
Vln.Ctrl_PCOS.list = list()
#celltype.x = "uNK 1"
for (celltype.x in plot.subtype) {
  
  print(paste("Plotting", celltype.x))
  
  # Get the seurat path with the celltype
  seurat.file = list.files(seurat.dir, full.names = T)[grep(celltype.x, list.files(seurat.dir))]
  
  # In case of error, HARDCODE
  if (length(seurat.file) == 0) {
    print("Resolve this manually")
    
    if (celltype.x == "SOX9+ LGR5+") {
      seurat.file = list.files(seurat.dir, full.names = T)[3]
    } else if (celltype.x == "SOX9+ proliferative") {
      seurat.file = list.files(seurat.dir, full.names = T)[4]
    }
    
  }
  
  # Load the seurat and extract the RNA normalised count matrix
  seurat.x = LoadH5Seurat(file = seurat.file)
  
  # Set the ident
  Idents(seurat.x) = seurat.ident
  
  # Reorder the ident to Group.order
  seurat.x$Group_Stage <- factor(x = seurat.x$Group_Stage, levels = Group.order)
  
  # Plot the violinplot of all groups
  Vln.x = VlnPlot(seurat.x, features = plot.genes, pt.size = 0, cols = Group.cols, 
                  group.by = "Group_Stage", assay = "RNA") + ggtitle(celltype.x)
  
  if (do.boxplot == TRUE) {
   Vln.x = Vln.x + geom_boxplot(width=0.3, fill="white", outlier.size = 0) 
  }
  
  # Attach the output
  Vln.all.list[[celltype.x]] = Vln.x
  
  # Save the Violinplot seperately
  ggsave2(filename = paste0(Output.dir, Project_name, celltype.x, "_Vlnplot_All.pdf"), plot = Vln.x, height = 4, width = 6)
  ggsave2(filename = paste0(Output.dir, Project_name, celltype.x, "_Vlnplot_All_small.pdf"), plot = Vln.x, height = 3, width = 4)
  
  # Plot individual genes if there are multiple genes to plot
  if (length(plot.genes) > 1) {
    
    for (gene.x in plot.genes) {
      
      # Plot the violinplot of all groups, single gene
      Vln.x = VlnPlot(seurat.x, features = gene.x, pt.size = 0, cols = Group.cols, 
                      group.by = "Group_Stage", assay = "RNA") + ggtitle(gene.x)
      Vln.x = Vln.x + geom_boxplot(width=0.2, fill="white", outlier.size = 0) 
      
      # Save the Violinplot seperately
      ggsave2(filename = paste0(Output.dir, Project_name, celltype.x, "_", gene.x, "_Vlnplot_All.pdf"), plot = Vln.x, height = 6, width = 10)
      
      # Attach the output
      Vln.all.gene.list[[paste0(celltype.x, "_", gene.x)]] = Vln.x
      
    } 
  }
  
  # Extra testing part
  #seurat.avg = AverageExpression(seurat.x, assays = "RNA")
  #avg.df = data.frame(seurat.avg$RNA)
  #avg.df$Gene = rownames(avg.df)
  
  #Dot.x = DotPlot(seurat.x, features = "PAEP", assay = "RNA", scale = FALSE)
  
  if (do.Ctrl_PCOS == TRUE) {
    # Plot the violinplot of only Ctrl and PCOS
    seurat.x = subset(seurat.x, idents = c("Control", "PCOS_W0"))
    Vln.y = VlnPlot(seurat.sub, features = plot.genes, pt.size = 0, cols = Group.cols, 
                    group.by = "Group_Stage", assay = "RNA")
    
    if (do.boxplot == TRUE) {
      Vln.y = Vln.y + geom_boxplot(width=0.3, fill="white") 
    }
    
    # Save the Violinplot seperately
    ggsave2(filename = paste0(Output.dir, Project_name, celltype.x, "_Vlnplot_Ctrl_PCOS.pdf"), plot = Vln.y, height = 6, width = 10)
    
    # Attach the output
    Vln.Ctrl_PCOS.list[[celltype.x]] = Vln.x
    
  }

}

# Do a gridplot of the plots
grid.Epithelium = plot_grid(plotlist=Vln.Ctrl_PCOS.list[1:3], ncol = 3)
grid.Stroma = plot_grid(plotlist=Vln.Ctrl_PCOS.list[4:6], ncol = 3)
grid.Epithelium.all = plot_grid(plotlist=Vln.all.list[1:3], ncol = 3)
grid.Stroma.all = plot_grid(plotlist=Vln.all.list[4:6], ncol = 3)
grid.Endo_Artery.up = plot_grid(plotlist = Vln.all.gene.list[1:6], ncol = 3, nrow = 2)
grid.Endo_Artery.down = plot_grid(plotlist = Vln.all.gene.list[7:8], ncol = 3, nrow = 2)
grid.Lymph = plot_grid(plotlist = Vln.all.gene.list[11], ncol = 3, nrow = 2)

grid.Endo_Artery.GO = plot_grid(plotlist = Vln.all.gene.list[1:4], ncol = 3, nrow = 2)
grid.Lymph.GO = plot_grid(plotlist = Vln.all.gene.list[c(9, 11:12)], ncol = 3, nrow = 2)

# Save the grid plots
ggsave2(filename = paste0(Output.dir, Project_name, "Vlnplot_Grid_Ctrl_PCOS_Epithelium.pdf"), plot = grid.Epithelium, height = 3, width = 8)
ggsave2(filename = paste0(Output.dir, Project_name, "Vlnplot_Grid_Ctrl_PCOS_Stromal.pdf"), plot = grid.Stroma, height = 3, width = 8)
ggsave2(filename = paste0(Output.dir, Project_name, "Vlnplot_Grid_All_Epithelium.pdf"), plot = grid.Epithelium.all, height = 3, width = 11)
ggsave2(filename = paste0(Output.dir, Project_name, "Vlnplot_Grid_All_Stromal.pdf"), plot = grid.Stroma.all, height = 3, width = 11)

ggsave2(filename = paste0(Output.dir, Project_name, "Vlnplot_Grid_Endo_Artery_Up.pdf"), plot = grid.Endo_Artery.up, height = 6, width = 11)
ggsave2(filename = paste0(Output.dir, Project_name, "Vlnplot_Grid_Endo_Artery_Down.pdf"), plot = grid.Endo_Artery.down, height = 6, width = 11)
ggsave2(filename = paste0(Output.dir, Project_name, "Vlnplot_Grid_Lymphatic.pdf"), plot = grid.Lymph, height = 6, width = 11)

ggsave2(filename = paste0(Output.dir, Project_name, "Vlnplot_Grid_Endo_Artery_GO.pdf"), plot = grid.Endo_Artery.GO, height = 6, width = 11)
ggsave2(filename = paste0(Output.dir, Project_name, "Vlnplot_Grid_Lymphatic_GO.pdf"), plot = grid.Lymph.GO, height = 6, width = 11)




