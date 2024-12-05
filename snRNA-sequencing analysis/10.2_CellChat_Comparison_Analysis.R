#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(CellChat)
library(openxlsx)
library(patchwork)
library(reticulate)
library(umap)
library(ComplexHeatmap)
library(wordcloud)
reticulate::use_python("/Users/gustaw.eriksson/anaconda3/bin/python", required=T)
#conda_install("numpy")


#Set work directory to directory where script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Generate or set the output directory
if (dir.exists(path = "Output/9_CellChat") == FALSE) {
  print("Output/9_CellChat")
  dir.create(path = "Output/9_CellChat", recursive = TRUE)
  Input.dir = "Output/9_CellChat/"
} else if (dir.exists(path = "Output/9_CellChat") == TRUE) {
  print("Directory exists")
  Input.dir = "Output/9_CellChat/"
} else {
  print("Error with output directory")
}

# Setting variables and parameters
Project_name = "Endo_All"
Group.cols = c("#A0A0A0", "#D098B2", "#95BFE1", "#65D46E")
Do.functional = FALSE
Do.structural = FALSE
Do.Merge_CellChat = TRUE
FC.cut = 0.5
padj.cut = 0.05
min.cell.exp = 0.25
gsub.pattern = ".*Endo_All_CellChat_Merged_object_(.*?).rds"
targets.x = c("Lumenal", "SOX9+ LGR5+", "SOX9+ LGR5-",
              "SOX9+ proliferative", "AR+", "Ciliated",
              "Stroma 1", "Stroma 2", 
              "Stroma proliferative", "Fibroblast", "uSMC",
              "T-cells CD4+", 
              "T-cells CD8+", "uNK 1", "uNK 2", "uNK 3", 
              "uM 1", "uM 2", "Tregs", "ILC3",                     
              "B-cells", "Mast cells", "DC1", "DC2", 
              "Migratory DC", "pDC",
              "Endothelial Vein", "Endothelial Artery", "Endothelial proliferative", "Mesenchymal",
              "Lymphatic")
endothelial.vec = c("Endothelial Vein", "Endothelial Artery", "Endothelial proliferative", "Mesenchymal",
                    "Lymphatic")
epithelium.vec = c("Lumenal", "SOX9+ LGR5+", "SOX9+ LGR5-",
                   "SOX9+ proliferative", "AR+", "Ciliated")
immune.vec = c("T-cells CD4+", 
               "T-cells CD8+", "uNK 1", "uNK 2", "uNK 3", 
               "uM 1", "uM 2", "Tregs", "ILC3",                     
               "B-cells", "Mast cells", "DC1", "DC2", 
               "Migratory DC", "pDC")
stroma.vec = c("Stroma 1", "Stroma 2", 
               "Stroma proliferative", "Fibroblast", "uSMC")
source.list = list(endothelial.vec, epithelium.vec, immune.vec, stroma.vec)

# Updated order and colors below:
subtype.order = c("Lumenal", "SOX9+ LGR5+", "SOX9+ LGR5-",              
                  "SOX9+ proliferative", "AR+", "Ciliated", # Epithelium
                  "Stroma 1", "Stroma 2", "Stroma proliferative", 
                  "Fibroblast", "uSMC", # Stroma
                  "Tregs", "T-cells CD8+", "T-cells CD4+",
                  "ILC3", "B-cells", "uNK 1", "uNK 2", "uNK 3",
                  "pDC", "Migratory DC", "DC1", "DC2", 
                  "uM 1", "uM 2", "Mast cells", # Immune
                  "Endothelial Vein", "Endothelial Artery", "Endothelial proliferative", 
                  "Mesenchymal", "Lymphatic") # Endothelial

subtype.col = c("#65C2A5", #Lumenal
                "#D3020D", #SOX9+ LGR5+
                "#2570B7", # SOX9+ LGR5-
                "#FAA0A1", # SOX9+ proliferative
                "#FF7F00", # AR+
                "#DD8D61", # Ciliated #DD8D61
                "#E64B35CC", # Stroma 1
                "#4DBBD599", # Stroma 2
                "#FAA0A6", # Stroma proliferative
                "#91D1C2CC", # Fibroblast
                "#F39B7FCC", # uSMC
                "salmon", #Treg
                "orangered", "darkred", # T-cells
                "violetred", # ILC3
                "sandybrown", # B-cells
                "dodgerblue", "skyblue", "navyblue", # uNKs
                "palegoldenrod", "rosybrown", # pDC and migratory DC
                "seagreen", "limegreen", # DC
                "plum", "maroon", # uMs
                "darkkhaki", #  mast cells
                "#8491B4CC", # Endothelial Vein
                "#B09C8599", # Endothelial Artery
                "#F39B7F99", # Endothelial proliferative
                "#00A08799", # Mesenchymal
                "#7E6148CC") # Lymphatic
names(subtype.col) = subtype.order

pathways.plot = c("COLLAGEN", "FN1", "LAMININ", "BMP", "IGF", 
                  "NRXN", "CNTN", "SLIT", "CADM", "SEMA3", "SPP1")

# List CellChat object to load
CC.list = list.files(path = Input.dir, pattern = paste0(Project_name, "_CellChat_object_.*.rds"), full.names = TRUE)
Group.names = sapply(CC.list, function(x) {gsub(".*_CellChat_object_(.*?).rds", "\\1", x)})

# Load the individual CellChat object per group
CC.list = lapply(CC.list, readRDS)
names(CC.list) = Group.names

# Editing order of vectors
CC.list = CC.list[c(1, 2, 4, 3)]
Group.names = Group.names[c(1, 2, 4, 3)]

# Reorder the CellChat object
CC.list = lapply(CC.list, function(x) updateClusterLabels(object = x, new.order = subtype.order))

# Load the combined and baseline CellChat object
if (Do.Merge_CellChat == TRUE) {
  # Merge the CellChat objects
  CC.merged.BL = mergeCellChat(CC.list[c(1,2)], add.names = c("Control", "PCOS_W0"))
  #CC.metformin = mergeCellChat(CC.list[c(2,4)], add.names = c("PCOS_W0", "PCOS_W16_Met"))
  #CC.lifestyle = mergeCellChat(CC.list[c(2,3)], add.names = c("PCOS_W0", "PCOS_W16_LS"))
  CC.merged.all <- mergeCellChat(CC.list, add.names = names(CC.list))
  CC.merged.list = list(CC.merged.all, CC.merged.BL)
} else if (Do.Merge_CellChat == FALSE) {
  
  # Load the already merged cellchats
  CC.merged.list = list.files(path = Input.dir, pattern = paste0(Project_name, ".*_CellChat_Merged_object_.*.rds"), full.names = TRUE)
  CC.merged.all = list.files(path = Input.dir, pattern = paste0(Project_name, "_CellChat_Merged_object_All.rds"), full.names = TRUE)
  CC.merged.BL = list.files(path = Input.dir, pattern = paste0(Project_name, "_CellChat_Merged_object_Control_PCOS.rds"), full.names = TRUE)
  CC.merged.list = CC.merged.list[CC.merged.list != CC.merged.all]
  
  CC.merged.all = readRDS(file = CC.merged.all)
  CC.merged.BL = readRDS(file = CC.merged.BL)
  
}

# Generate or set the output directory
if (dir.exists(path = "Output/9_CellChat/Comparison_analysis") == FALSE) {
  print("Output/9_CellChat/Comparison_analysis")
  dir.create(path = "Output/9_CellChat/Comparison_analysis", recursive = TRUE)
  Output.dir = "Output/9_CellChat/Comparison_analysis/"
} else if (dir.exists(path = "Output/9_CellChat/Comparison_analysis") == TRUE) {
  print("Directory exists")
  Output.dir = "Output/9_CellChat/Comparison_analysis/"
} else {
  print("Error with output directory")
}

# Run analysis on all samples in mergead all dataset
#Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(CC.merged.all, show.legend = F, group = c(1:length(CC.merged.all@var.features)), color.use = Group.cols[c(1:4)])
gg2 <- compareInteractions(CC.merged.all, show.legend = F, group = c(1:length(CC.merged.all@var.features)), measure = "weight", color.use = Group.cols[c(1:4)])
gg.plot = gg1 + gg2
ggsave(filename = paste0(Output.dir, "Compare_interaction_All_barplot.pdf"), plot = gg.plot, width = 10, height = 7)

# Run netAnalysis_computeCentrality on all cellchat objects
#CC.list <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
CC.list = lapply(CC.list, netAnalysis_computeCentrality)

# Function to compare the major sources and targets in 2D space
Celltype_2D_plot <- function(CC.plot.list, x.name, axis.lim = c(0, 0.45)) {
  num.link <- sapply(CC.plot.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(CC.plot.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(CC.plot.list[[i]], title = names(CC.plot.list)[i], weight.MinMax = weight.MinMax, color.use = subtype.col) + 
      xlim(axis.lim) + ylim(axis.lim)
  }
  gg.plot = patchwork::wrap_plots(plots = gg)
  ggsave(filename = paste0(Output.dir, "Compare_Source_Target_2D_", x.name, "_plot.pdf"), plot = gg.plot, width = 8, height = 8)
  return(gg.plot)
}

All_Celltype_2D = Celltype_2D_plot(CC.plot.list = CC.list, x.name = "All")
BL_Celltype_2D = Celltype_2D_plot(CC.plot.list = CC.list[c(1,2)], x.name = "BL")
All_Celltype_2D
BL_Celltype_2D

# NetAnalysis scatter plot on all cell types
Celltype_Comp_2D_plot <- function(celltype.x, comp.x = 1, comp.y = 2, comp.name = "Ctrl_PCOS", axis.lim = c(0, 0.05)) {
  
  if (is.null(axis.lim) == TRUE) {
    gg.plot <- netAnalysis_signalingChanges_scatter(CC.list, idents.use = celltype.x, color.use = subtype.col, comparison = c(comp.x, comp.y))
  } else if (is.null(axis.lim) == FALSE) {
    gg.plot <- netAnalysis_signalingChanges_scatter(CC.list, idents.use = celltype.x, color.use = subtype.col, comparison = c(comp.x, comp.y)) + 
      xlim(axis.lim) + ylim(axis.lim)
  }
  ggsave(filename = paste0(Output.dir, "Compare_Source_Target_2D_", celltype.x, "_", comp.name, "_plot.pdf"), plot = gg.plot, width = 8, height = 8)
  return(gg1)
}

Celltype_comp_2D_plots = lapply(subtype.order, function(celltype) Celltype_Comp_2D_plot(celltype.x = celltype, axis.lim = NULL))
Celltype_comp_2D_plots = lapply(subtype.order, function(celltype) Celltype_Comp_2D_plot(celltype.x = celltype, axis.lim = NULL, comp.x = 2, comp.y = 3, comp.name = "PCOS_Metformin"))
Celltype_comp_2D_plots = lapply(subtype.order, function(celltype) Celltype_Comp_2D_plot(celltype.x = celltype, axis.lim = NULL, comp.x = 2, comp.y = 4, comp.name = "PCOS_Lifestyle"))

# Computing netsimilarity takes time. As we have the same cell types between the datasets, functional analysis is recommended
if (Do.functional == TRUE) {
  
  #Compute functional similarity only Control vs. PCOS
  # Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
  CC.merged.BL <- computeNetSimilarityPairwise(CC.merged.BL, type = "functional")
  CC.merged.BL <- netEmbedding(CC.merged.BL, type = "functional")
  CC.merged.BL <- netClustering(CC.merged.BL, type = "functional")
  
  # Visualization in 2D-space
  gg.plot = netVisual_embeddingPairwise(CC.merged.BL, type = "functional", label.size = 3.5)
  ggsave(filename = paste0(Output.dir, "Functional_Similarity_2D_BL_plot.pdf"), plot = gg.plot, width = 10, height = 7)
  
  netVisual_embeddingZoomIn(CC.merged.BL, type = "functional", nCol = 2)
  
  gg.plot = netVisual_embeddingPairwise(CC.merged.BL, type = "functional", label.size = 3.5, 
                                        do.label = TRUE, pathway.labeled = "JAM")
  gg.plot
  ggsave(filename = paste0(Output.dir, "Functional_Similarity_2D_BL_Selectedplot.pdf"), plot = gg.plot, width = 10, height = 7)
  
  netVisual_embedding(CC.merged.BL, type = "functional")
  
  gg.plot = rankSimilarity(CC.merged.BL, type = "functional", bar.w = 0.9) + 
    geom_bar(stat = "identity", fill = "#669999", ) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme_cowplot()
  gg.plot
  ggsave(filename = paste0(Output.dir, "Functional_Similarity_Barplot_BL_plot.pdf"), plot = gg.plot, width = 5, height = 10)
  
  # Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
  CC.merged.all <- computeNetSimilarityPairwise(CC.merged.all, type = "functional")
  CC.merged.all <- netEmbedding(CC.merged.all, type = "functional")
  CC.merged.all <- netClustering(CC.merged.all, type = "functional")
  
  # Visualization in 2D-space
  gg.plot = netVisual_embeddingPairwise(CC.merged.all, type = "functional", label.size = 3.5)
  ggsave(filename = paste0(Output.dir, "Functional_Similarity_2D_All_plot.pdf"), plot = gg.plot, width = 10, height = 7)
  gg.plot = rankSimilarity(CC.merged.all, type = "functional", bar.w = 0.9) + 
    geom_bar(stat = "identity", fill = "#669999", ) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme_cowplot()
  gg.plot
  ggsave(filename = paste0(Output.dir, "Functional_Similarity_Barplot_All_plot.pdf"), plot = gg.plot, width = 5, height = 10)
}

if (Do.structural == TRUE) {
  CC.merged.all <- computeNetSimilarityPairwise(CC.merged.all, type = "structural")
  CC.merged.all <- netEmbedding(CC.merged.all, type = "structural")
  CC.merged.all <- netClustering(CC.merged.all, type = "structural")
  # Visualization in 2D-space
  gg.plot = netVisual_embeddingPairwise(CC.merged.all, type = "structural", label.size = 3.5)
  ggsave(filename = paste0(Output.dir, "Structural_Similarity_2D_All_plot.pdf"), plot = gg.plot, width = 10, height = 7)
}

CC.x = CC.merged.list[1]
source.x = epithelium.vec
target.x = targets.x

RankNet_plot <- function(CC.x, source.x, target.x, name.x) {
  
  # Read the CellChat object
  CC.x = readRDS(CC.x)
  
  # Create output list
  Output.list = list()
  
  # Generate the RankNetand save output
  gg1 = rankNet(CC.x, sources.use = source.x, targets.use = target.x, stacked = TRUE, do.stat = TRUE, return.data = TRUE, color.use = Group.cols[c(1,2)])
  #gg1$gg.obj
  #ggsave(plot = gg1$gg.obj, filename = paste0(Output.dir, paste0("RankNet_barplot_", name.x, ".pdf")), plot = gg.plot, width = 10, height = 7)
  Output.list$Source = gg1
  
  # Generate the RankNetand save output
  gg1 = rankNet(CC.x, sources.use = target.x, targets.use = source.x, stacked = TRUE, do.stat = TRUE, return.data = TRUE, color.use = Group.cols[c(1,2)])
  #gg1$gg.obj
  #ggsave(plot = gg1$gg.obj, filename = paste0(Output.dir, paste0("RankNet_barplot_", name.x, ".pdf")), plot = gg.plot, width = 10, height = 7)
  Output.list$TargetvsSource = gg1
  
  # Loop over each celltype and test its interacation to target
  for (celltype.x in source.x) {
    
    gg1 = rankNet(CC.x, sources.use = celltype.x, targets.use = target.x, stacked = TRUE, do.stat = TRUE, return.data = TRUE, color.use = Group.cols[c(1,2)])
    #gg1$gg.obj
    #ggsave(plot = gg1$gg.obj, filename = paste0(Output.dir, paste0("RankNet_barplot_", name.x, "_", celltype.x, ".pdf")), plot = gg.plot, width = 10, height = 7)
    Output.list[[celltype.x]] = gg1
    
  }
  
  return(Output.list)
  
} 

RankNet_epithelium = RankNet_plot(CC.x = CC.merged.list[1], source.x = epithelium.vec, target.x = targets.x, name.x = "Ctrl_PCOS_Epithelium")
RankNet_stroma = RankNet_plot(CC.x = CC.merged.list[1], source.x = stroma.vec, target.x = targets.x, name.x = "Ctrl_PCOS_Stroma")
RankNet_immune = RankNet_plot(CC.x = CC.merged.list[1], source.x = immune.vec, target.x = targets.x, name.x = "Ctrl_PCOS_Immune")
rankNet_endothelium = RankNet_plot(CC.x = CC.merged.list[1], source.x = endothelial.vec, target.x = targets.x, name.x = "Ctrl_PCOS_Endothelium")
rankNet_all = RankNet_plot(CC.x = CC.merged.list[1], source.x = targets.x, target.x = targets.x, name.x = "Ctrl_PCOS_All")

# Function to plot proportion of variable in dataframe
Prop_barplot <- function(x, x.name = "Up", x.col = c("Ligand", "Percantage"), CC.df, plot.color = "black") {
  
  # Generate the proportion table
  CC.prop = table(CC.df[[x]])
  CC.prop = as.data.frame(prop.table(CC.prop)*100)
  CC.prop = CC.prop[order(CC.prop[,2], decreasing = FALSE), ]
  colnames(CC.prop) = x.col
  
  # Generate a barplot of the proportions
  gg.plot = ggplot(CC.prop, aes(x = factor(CC.prop[,1], levels = CC.prop[,1]), y = CC.prop[,2])) +
    geom_bar(stat = "identity", fill = plot.color) +
    labs(x = x.col[1], y = x.col[2]) + 
    theme_cowplot() + theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    coord_flip()
  
  ggsave(filename = paste0(Output.dir, "Prop_Barplot_", x.name, "_", x, ".pdf"), plot = gg.plot, width = 10, height = 7)
  
  return(gg.plot)
  
}

Run_Comparison_analysis <- function(CC.comp, CC.i, CC.j, pos.dataset = "PCOS_W0", CC.type = "Control_PCOS",
                                    read.RDS = FALSE) {
  
  # Print what is being analysed
  print(paste("Running CellChat on", CC.type))
  
  # Load the CellChat merged object
  if (read.RDS == TRUE) {
    CC.x = readRDS(CC.comp) 
  } else if (read.RDS == FALSE) {
    CC.x = CC.comp
    CC.comp = NULL
  }
  
  # Plotting number of interactions and interaction strength in selected comparison
  gg1 <- compareInteractions(CC.x, show.legend = F, group = c(1:2), color.use = Group.cols[c(CC.i, CC.j)])
  gg2 <- compareInteractions(CC.x, show.legend = F, group = c(1:2), measure = "weight", color.use = Group.cols[c(CC.i, CC.j)])
  plot.all = gg1 + gg2
  ggsave(filename = paste0(Output.dir, "CompareInteraction_", CC.type, "_plot.pdf"), plot = plot.all)

  # Generate heatmaps of number of interaction and strength
  gg1 <- netVisual_heatmap(CC.x, color.use = subtype.col)
  gg2 <- netVisual_heatmap(CC.x, measure = "weight", color.use = subtype.col)
  pdf(file=paste0(Output.dir, "CompareInteractions_heatmap_Ligand_Receptor_", CC.type, ".pdf"), width = 14, height = 9)
  print(gg1 + gg2)
  dev.off()
  
  # Identify and visualize the conserved and context-specific signaling pathways
  ## Compare the overall information flow of each signaling pathway
  gg1 <- rankNet(CC.x, measure = "count", mode = "comparison", stacked = T, do.stat = TRUE, color.use = Group.cols[c(CC.i, CC.j)], return.data = F)
  gg2 <- rankNet(CC.x, measure = "weight", mode = "comparison", stacked = T, do.stat = TRUE, color.use = Group.cols[c(CC.i, CC.j)], return.data = F)
  pdf(file=paste0(Output.dir, "RankNet_Information_flow_barplot_", CC.type, ".pdf"), height = 9)
  print(gg1 + gg2)
  dev.off()

  ######### Generate a list of the individual Cell Chat object called "object list" ######
  # Compare outgoing (and all or incoming) signaling associated with each cell population
  # combining all the identified signaling pathways from different datasets 
  p.union <- union(CC.list[[CC.i]]@netP$pathways, CC.list[[CC.j]]@netP$pathways)
  ht1 = netAnalysis_signalingRole_heatmap(CC.list[[CC.i]], pattern = "outgoing", signaling = p.union, title = names(CC.list)[CC.i], width = 8, height = 20, color.use = subtype.col)
  ht2 = netAnalysis_signalingRole_heatmap(CC.list[[CC.j]], pattern = "outgoing", signaling = p.union, title = names(CC.list)[CC.j], width = 8, height = 20, color.use = subtype.col)
  pdf(file=paste0(Output.dir, "Outgoing_signal_heatmap_", CC.type, ".pdf"), height = 12, width = 10)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  ht1 = netAnalysis_signalingRole_heatmap(CC.list[[CC.i]], pattern = "incoming", signaling = p.union, title = names(CC.list)[CC.i], width = 8, height = 20, color.use = subtype.col)
  ht2 = netAnalysis_signalingRole_heatmap(CC.list[[CC.j]], pattern = "incoming", signaling = p.union, title = names(CC.list)[CC.j], width = 8, height = 20, color.use = subtype.col)
  pdf(file=paste0(Output.dir, "Incoming_signal_heatmap_", CC.type, ".pdf"), height = 12, width = 10)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  ht1 = netAnalysis_signalingRole_heatmap(CC.list[[CC.i]], pattern = "all", signaling = p.union, title = names(CC.list)[CC.i], width = 8, height = 20, color.use = subtype.col, cluster.cols = T)
  ht2 = netAnalysis_signalingRole_heatmap(CC.list[[CC.j]], pattern = "all", signaling = p.union, title = names(CC.list)[CC.j], width = 8, height = 20, color.use = subtype.col, cluster.cols = T)
  pdf(file=paste0(Output.dir, "All_signal_heatmap_", CC.type, ".pdf"), height = 12, width = 10)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
  #######
  
  # Generate heatmpas with selected pathways
  p.union.select = p.union[p.union %in% pathways.plot]
  ht1 = netAnalysis_signalingRole_heatmap(CC.list[[CC.i]], pattern = "outgoing", signaling = p.union.select, title = names(CC.list)[CC.i], width = 9, height = 3, color.use = subtype.col)
  ht2 = netAnalysis_signalingRole_heatmap(CC.list[[CC.j]], pattern = "outgoing", signaling = p.union.select, title = names(CC.list)[CC.j], width = 9, height = 3, color.use = subtype.col)
  pdf(file=paste0(Output.dir, "Outgoing_signal_heatmap_", CC.type, "_Selected.pdf"), height = 6, width = 10)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  ht1 = netAnalysis_signalingRole_heatmap(CC.list[[CC.i]], pattern = "incoming", signaling = p.union.select, title = names(CC.list)[CC.i], width = 9, height = 3, color.use = subtype.col)
  ht2 = netAnalysis_signalingRole_heatmap(CC.list[[CC.j]], pattern = "incoming", signaling = p.union.select, title = names(CC.list)[CC.j], width = 9, height = 3, color.use = subtype.col)
  pdf(file=paste0(Output.dir, "Incoming_signal_heatmap_", CC.type, "_Selected.pdf"), height = 12, width = 10)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  ht1 = netAnalysis_signalingRole_heatmap(CC.list[[CC.i]], pattern = "all", signaling = p.union.select, title = names(CC.list)[CC.i], width = 9, height = 3, color.use = subtype.col)
  ht2 = netAnalysis_signalingRole_heatmap(CC.list[[CC.j]], pattern = "all", signaling = p.union.select, title = names(CC.list)[CC.j], width = 9, height = 3, color.use = subtype.col)
  pdf(file=paste0(Output.dir, "All_signal_heatmap_", CC.type, "_Selected.pdf"), height = 12, width = 10)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  #Compare the major sources and targets in 2D space
  Comp_2D_plot = Celltype_2D_plot(CC.plot.list = CC.list[c(CC.i, CC.j)], x.name = CC.type, axis.lim = c(0, 0.5))
  
  # Do Scatterplots of differentually signalling in each celltype
  celltype.x = subtype.order[2]
  for (celltype.x in subtype.order) {
    
    gg.plot <- netAnalysis_signalingChanges_scatter(CC.x, idents.use = celltype.x, color.use = subtype.col)
    ggsave(filename = paste0(Output.dir, "SignalingChange_Scatter_", CC.type, "_", celltype.x, "_plot.pdf"), plot = gg.plot, width = 5, height = 3)
    
    
  }
  
  # Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
  CC.x <- computeNetSimilarityPairwise(CC.x, type = "functional")
  CC.x <- netEmbedding(CC.x, type = "functional")
  CC.x <- netClustering(CC.x, type = "functional")
  # Visualization in 2D-space
  gg.plot = netVisual_embeddingPairwise(CC.x, type = "functional", label.size = 3.5)
  ggsave(filename = paste0(Output.dir, "Functional_Similarity_2D_", CC.type, "_plot.pdf"), plot = gg.plot, width = 10, height = 7)
  gg.plot = rankSimilarity(CC.x, type = "functional", bar.w = 0.9) + 
    geom_bar(stat = "identity", fill = "#669999", ) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme_cowplot()
  gg.plot
  ggsave(filename = paste0(Output.dir, "Functional_Similarity_Barplot_", CC.type, "_plot.pdf"), plot = gg.plot, width = 5, height = 10)
  
  CC.x <- computeNetSimilarityPairwise(CC.x, type = "structural")
  CC.x <- netEmbedding(CC.x, type = "structural")
  CC.x <- netClustering(CC.x, type = "structural")
  # Visualization in 2D-space
  gg.plot = netVisual_embeddingPairwise(CC.x, type = "structural", label.size = 3.5)
  ggsave(filename = paste0(Output.dir, "Structural_Similarity_2D_", CC.type, "_plot.pdf"), plot = gg.plot, width = 10, height = 7)
  gg.plot = rankSimilarity(CC.x, type = "structural", bar.w = 0.9) + 
    geom_bar(stat = "identity", fill = "#669999", ) + 
    scale_y_continuous(expand = c(0, 0)) +
    theme_cowplot()
  gg.plot
  ggsave(filename = paste0(Output.dir, "Structural_Similarity_Barplot_", CC.type, "_plot.pdf"), plot = gg.plot, width = 5, height = 10)
  
  # Circle plot of differential number of interactions and strength
  pdf(file=paste0(Output.dir, "netVisual_", CC.type, "_Difference_Interaction_Number", ".pdf"))
  netVisual_diffInteraction(CC.x, weight.scale = T, color.use = subtype.col)
  dev.off()
  
  pdf(file=paste0(Output.dir, "netVisual_", CC.type, "_Difference_Interaction_Strenght", ".pdf"))
  netVisual_diffInteraction(CC.x, weight.scale = T, measure = "weight", color.use = subtype.col)
  dev.off()
  
  # Netvisual if 
  gg1 <- netVisual_heatmap(CC.x, color.use = subtype.col, height = 8, width = 4)
  gg2 <- netVisual_heatmap(CC.x, measure = "weight", color.use = subtype.col, height = 8, width = 4)
  pdf(file=paste0(Output.dir, "CompareInteractions_heatmap_Ligand_Receptor_", CC.type, ".pdf"), width = 14, height = 10)
  print(gg1 + gg2)
  dev.off()
  
  # Compare cell-cell communication of specific pathways
  for (pathway.x in pathways.plot) {
    
    # Do circle plot
    print(pathway.x)
    weight.max <- getMaxWeight(CC.list[c(CC.i, CC.j)], slot.name = c("netP"), attribute = pathway.x) # control the edge weights across different datasets
    pdf(file=paste0(Output.dir, "netVisual_aggregate_", CC.type, "_CirclePlot_", names(CC.list[CC.i]), "_",  pathway.x, ".pdf"))
    netVisual_aggregate(CC.list[[CC.i]], signaling = pathway.x, layout = "circle", edge.weight.max = weight.max[1], 
                        edge.width.max = 10, signaling.name = paste(pathway.x, names(CC.list)[CC.i]), color.use = subtype.col)
    dev.off()
    pdf(file=paste0(Output.dir, "netVisual_aggregate_", CC.type, "_CirclePlot_", names(CC.list[CC.j]), "_",  pathway.x, ".pdf"))
    netVisual_aggregate(CC.list[[CC.j]], signaling = pathway.x, layout = "circle", 
                        edge.weight.max = weight.max[1], edge.width.max = 10, 
                        signaling.name = paste(pathway.x, names(CC.list)[CC.j]), , color.use = subtype.col)
    dev.off()
    
    # Do chard plot
    pdf(file=paste0(Output.dir, "netVisual_aggregate_", CC.type, "_ChordPlot_", names(CC.list[CC.i]), "_",  pathway.x, ".pdf"))
    netVisual_aggregate(CC.list[[CC.i]], signaling = pathway.x, layout = "chord", 
                        signaling.name = paste(pathway.x, names(CC.list)[CC.i]), color.use = subtype.col)
    dev.off()
    pdf(file=paste0(Output.dir, "netVisual_aggregate_", CC.type, "_ChordPlot_", names(CC.list[CC.j]), "_",  pathway.x, ".pdf"))
    netVisual_aggregate(CC.list[[CC.j]], signaling = pathway.x, layout = "chord", 
                        signaling.name = paste(pathway.x, names(CC.list)[CC.j]), , color.use = subtype.col)
    dev.off()
    

    #weight.max <- getMaxWeight(CC.list, slot.name = c("netP"), attribute = pathway.x) # control the edge weights across different datasets
    #par(mfrow = c(1,4), xpd=TRUE)
    #for (i in 1:length(CC.list)) {
    #  netVisual_aggregate(CC.list[[i]], signaling = pathway.x, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(CC.list)[i]))
    #}
    
  }
      
  # Identify dysfunctional signaling by using differential expression analysis
  ## define a char name used for storing the results of differential expression analysis
  features.name = pos.dataset
  
  # Map the results of differential expression analysis onto the inferred cell-cell communications 
  # to easily manage/subset the ligand-receptor pairs of interest
  CC.net <- netMappingDEG(CC.x, features.name = features.name)
  CC.up.ligand <- subsetCommunication(CC.x, net = CC.net, datasets = features.name, thresh = 0.05, ligand.logFC = 0.5)
  CC.up.receptor <- subsetCommunication(CC.x, net = CC.net, datasets = features.name, thresh = 0.05, receptor.logFC = 0.5)
  CC.down.ligand <- subsetCommunication(CC.x, net = CC.net, datasets = features.name, thresh = 0.05, ligand.logFC = -0.5)
  CC.down.receptor <- subsetCommunication(CC.x, net = CC.net, datasets = features.name, thresh = 0.05, receptor.logFC = -0.5)
  
  # Combine the vectors

  CC.DE.analysis = list(CC.net, CC.up.ligand, CC.up.receptor, CC.down.ligand, CC.down.receptor)
  names(CC.DE.analysis) = c("All_DE", "Up_ligand", "Up_receptor", "Down_ligand", "Down_receptor")
  
  # Save each list object as an Excel file
  lapply(names(CC.DE.analysis), function(name) {
    write.xlsx(CC.DE.analysis[[name]], file = paste0(Output.dir, name, "_table_", CC.type, ".xlsx"))
  })
  
  return(CC.DE.analysis)

}

# Run the comparison analysis, outputting differentially expressed ligands and receptors
Comp.Ctrl_vs_PCOS = Run_Comparison_analysis(CC.comp = CC.merged.list[[2]], CC.i = 1, CC.j = 2, pos.dataset = "PCOS_W0", CC.type = "Control_PCOS")
Comp.Ctrl_vs_PCOS = Run_Comparison_analysis(CC.comp = CC.merged.list[2], CC.i = 1, CC.j = 2, pos.dataset = "PCOS_W0")
Comp.PCOS_vs_LS = Run_Comparison_analysis(CC.comp = CC.merged.list[2], CC.i = 2, CC.j = 4, pos.dataset = "PCOS_W16_LS")
Comp.PCOS_vs_Met = Run_Comparison_analysis(CC.comp = CC.merged.list[3], CC.i = 2, CC.j = 3, pos.dataset = "PCOS_W0")

# Save the output
saveRDS(Comp.Ctrl_vs_PCOS, file = paste0(Output.dir, "Comparison_analysis_output_Ctrl_vs_PCOS.rds"))
saveRDS(Comp.PCOS_vs_LS, file = paste0(Output.dir, "Comparison_analysis_output_PCOS_vs_LS.rds"))
saveRDS(Comp.PCOS_vs_Met, file = paste0(Output.dir, "Comparison_analysis_output_PCOS_vs_Met.rds"))
