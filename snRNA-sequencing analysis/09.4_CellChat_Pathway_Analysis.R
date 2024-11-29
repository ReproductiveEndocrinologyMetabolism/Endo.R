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
library(RColorBrewer)
library(NMF)
library(ggalluvial)

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
Group.cols = c("#A0A0A0", "#D098B2", "#65D46E", "#95BFE1")
gsub.pattern = ".*Endo_All_CellChat_object_(.*?).rds"
targets.x = c("Endothelial Vein", "Endothelial Artery", "Endothelial proliferative", "Mesenchymal",
              "Lymphatic", "Lumenal", "SOX9+ LGR5+", "SOX9+ LGR5-",
              "SOX9+ proliferative", "AR+", "Ciliated", "T-cells CD4+", 
              "T-cells CD8+", "uNK 1", "uNK 2", "uNK 3", 
              "uM 1", "uM 2", "Tregs", "ILC3",                     
              "B-cells", "Mast cells", "DC1", "DC2", 
              "Migratory DC", "pDC", "Stroma 1", "Stroma 2", 
              "Stroma proliferative", "Fibroblast", "uSMC")

# Set the number of Kmers from running selectK
#CC.selectK.out = c(5, 3, 4, 2) # Maybe 4 on metformin # Order of the sample groups
#CC.selectK.in = c(3, 5, 2, 3) # Order of the sample groups
#run.selectK = FALSE

# Updated Kmers
CC.selectK.out = c(4, 3, 3, 2) # Maybe 4 on metformin # Order of the sample groups
CC.selectK.in = c(3, 3, 2, 3) # Order of the sample groups
run.selectK = FALSE

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
names(subtype.col_alt) = subtype.order_alt

#pathways.plot = c("COLLAGEN", "FN1", "LAMININ", "MK", "THBS")
pathways.plot = c("COLLAGEN", "FN1", "LAMININ", "BMP", "IGF", 
                  "NRXN", "CNTN", "SLIT", "CADM", "SEMA3", "SPP1")

# List CellChat object to load
CC.list = list.files(path = Input.dir, pattern = paste0(Project_name, ".*_CellChat_object_.*.rds"), full.names = TRUE)
Group.names = sapply(CC.list, function(x) {gsub(".*_CellChat_object_(.*?).rds", "\\1", x)})

# Load the individual CellChat object per group
CC.list = lapply(CC.list, readRDS)
names(CC.list) = Group.names

# Editing order of vectors
#names(CC.list) = names(CC.list)[c(1, 2, 4, 3)]
#Group.names = Group.names[c(1, 2, 4, 3)]
#names(Group.cols) = Group.names

# Generate or set the output directory
if (dir.exists(path = "Output/9_CellChat/Targeted_analysis") == FALSE) {
  print("Output/9_CellChat/Targeted_analysis")
  dir.create(path = "Output/9_CellChat/Targeted_analysis", recursive = TRUE)
  Output.dir = "Output/9_CellChat/Targeted_analysis/"
} else if (dir.exists(path = "Output/9_CellChat/Targeted_analysis") == TRUE) {
  print("Directory exists")
  Output.dir = "Output/9_CellChat/Targeted_analysis/"
} else {
  print("Error with output directory")
}

CC.list = lapply(CC.list, function(x) updateClusterLabels(object = x, new.order = subtype.order_alt))

# Loop over the CC objects
for (CC.i in 1:length(CC.list)) {
  
  # Extract the name of the CellCHat object
  CC.name = names(CC.list)[CC.i]
  print(CC.name)
  
  # Extract the current CC object
  CC.x = CC.list[[CC.i]]
  
  # Reorder the object
  #updateClusterLabels(CC.x, new.order = subtype.order_alt)
  #netAnalysis_computeCentrality(CC.x, slot.name = "netP")
  
  # Generate a dataframe of the results
  CC.df = subsetCommunication(CC.x)
  openxlsx::write.xlsx(x = CC.df, file = paste0(Output.dir, "Communication_table_", CC.name, ".xlsx"))
  groupSize <- as.numeric(table(CC.x@idents))
  gg1 = netVisual_circle(CC.x@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",
                        color.use = subtype.col_alt)
  pdf(file=paste0(Output.dir, "netVisual_", CC.name, "_Interaction_Strenght", ".pdf"))
  print(gg1)
  dev.off()
  
  # Extract most common pathways
  CC.kegg = table(CC.df$evidence)
  CC.kegg.prop = as.data.frame(prop.table(CC.kegg) * 100)
  CC.kegg = as.data.frame(CC.kegg)
  CC.kegg = as.data.frame(table(CC.df$evidence))
  colnames(CC.kegg) <- c("Category", "Count")
  colnames(CC.kegg.prop) <- c("Category", "Percentage")
  CC.kegg = CC.kegg[order(CC.kegg$Count, decreasing = TRUE), ]
  CC.kegg.10 = head(CC.kegg[order(-CC.kegg$Count), ], 10)
  CC.kegg.prop = CC.kegg.prop[order(CC.kegg.prop$Percentage, decreasing = TRUE), ]
  CC.kegg.prop.10 = head(CC.kegg.prop[order(-CC.kegg.prop$Percentage), ], 10)
  
  gg1 = ggplot(CC.kegg.10, aes(x = factor(Category, levels = CC.kegg.10$Category), y = Count)) +
    geom_bar(stat = "identity", fill = Group.cols[CC.i]) +
    labs(title = paste0("KEGG in ", CC.name), x = "Categories", y = "Count") + 
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  gg2 = ggplot(CC.kegg.prop.10, aes(x = factor(Category, levels = CC.kegg.prop.10$Category), y = Percentage)) +
    geom_bar(stat = "identity", fill = Group.cols[CC.i]) +
    labs(title = paste0("KEGG in ", CC.name), x = "Categories", y = "Proportion %") + 
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave2(filename = paste0(Output.dir, "Prop_Quant_pathway_barplot_", CC.name, ".pdf"), plot = gg1 + gg2)
  
  # Extract most common cell source
  CC.source = table(CC.df$source)
  CC.source.prop = as.data.frame(prop.table(CC.source) * 100)
  CC.source = as.data.frame(CC.source)
  CC.source = as.data.frame(table(CC.df$source))
  colnames(CC.source) <- c("Category", "Count")
  colnames(CC.source.prop) <- c("Category", "Percentage")
  CC.source = CC.source[order(CC.source$Count, decreasing = TRUE), ]
  CC.source.10 = head(CC.source[order(-CC.source$Count), ], 10)
  CC.source.prop = CC.source.prop[order(CC.source.prop$Percentage, decreasing = TRUE), ]
  CC.source.prop.10 = head(CC.source.prop[order(-CC.source.prop$Percentage), ], 10)
  
  gg1 = ggplot(CC.source.prop, aes(x = factor(Category, levels = CC.source.prop$Category), y = Percentage)) +
    geom_bar(stat = "identity", fill = Group.cols[CC.i]) +
    labs(title = paste0("Bar Plot source in ", CC.name), x = "Categories", y = "Proportion %") + 
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave2(filename = paste0(Output.dir, "Prop_celltype_barplot_", CC.name, ".pdf"), plot = gg1)
  
  # Extract most common pathway name
  CC.pathway_name = table(CC.df$pathway_name)
  CC.pathway_name.prop = as.data.frame(prop.table(CC.pathway_name) * 100)
  CC.pathway_name = as.data.frame(CC.pathway_name)
  CC.pathway_name = as.data.frame(table(CC.df$pathway_name))
  colnames(CC.pathway_name) <- c("Category", "Count")
  colnames(CC.pathway_name.prop) <- c("Category", "Percentage")
  CC.pathway_name = CC.pathway_name[order(CC.pathway_name$Count, decreasing = TRUE), ]
  CC.pathway_name.10 = head(CC.pathway_name[order(-CC.pathway_name$Count), ], 10)
  CC.pathway_name.20 = head(CC.pathway_name[order(-CC.pathway_name$Count), ], 20)
  CC.pathway_name.30 = head(CC.pathway_name[order(-CC.pathway_name$Count), ], 30)
  CC.pathway_name.40 = head(CC.pathway_name[order(-CC.pathway_name$Count), ], 40)
  CC.pathway_name.prop = CC.pathway_name.prop[order(CC.pathway_name.prop$Percentage, decreasing = TRUE), ]
  CC.pathway_name.prop.10 = head(CC.pathway_name.prop[order(-CC.pathway_name.prop$Percentage), ], 10)
  CC.pathway_name.prop.20 = head(CC.pathway_name.prop[order(-CC.pathway_name.prop$Percentage), ], 20)
  CC.pathway_name.prop.30 = head(CC.pathway_name.prop[order(-CC.pathway_name.prop$Percentage), ], 30)
  CC.pathway_name.prop.40 = head(CC.pathway_name.prop[order(-CC.pathway_name.prop$Percentage), ], 40)
  
  gg1 = ggplot(CC.pathway_name.prop, aes(x = factor(Category, levels = CC.pathway_name.prop$Category), y = Percentage)) +
    geom_bar(stat = "identity", fill = Group.cols[CC.i]) +
    labs(title = paste0("Bar Plot pathway in ", CC.name), x = "Categories", y = "Proportion %") + 
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ggsave2(filename = paste0(Output.dir, "Prop_pathway_barplot_", CC.name, ".pdf"), plot = gg1, width = 10)
  
  # Compute the network centrality scores
  CC.x = netAnalysis_computeCentrality(CC.x, slot.name = "netP")
  gg1 <- netAnalysis_signalingRole_scatter(CC.x, color.use = subtype.col)
  pdf(file=paste0(Output.dir, "Scatterplot_CentralityScore_", CC.name, "_All_Pathways.pdf"), width = 9, height = 7)
  print(gg1)
  dev.off()
  
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  gg1 <- netAnalysis_signalingRole_heatmap(CC.x, pattern = "outgoing", height = 18, 
                                           color.use = subtype.col)
  gg2 <- netAnalysis_signalingRole_heatmap(CC.x, pattern = "incoming", height = 18, 
                                           color.use = subtype.col)
  gg1 = gg1 + gg2
  pdf(file=paste0(Output.dir, "HeatmaP_CentralityScore_InOutSignaling_", CC.name, "_All_pathways.pdf"), width = 12, height = 18)
  print(gg1)
  dev.off()
  
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  gg1 <- netAnalysis_signalingRole_heatmap(CC.x, pattern = "outgoing", height = 18,
                                           color.use = subtype.col, signaling = CC.pathway_name.30$Category)
  gg2 <- netAnalysis_signalingRole_heatmap(CC.x, pattern = "incoming", height = 18,
                                           color.use = subtype.col, signaling = CC.pathway_name.prop.30$Category)
  gg1 = gg1 + gg2
  pdf(file=paste0(Output.dir, "HeatmaP_CentralityScore_InOutSignaling_", CC.name, "_Top30_pathways.pdf"), width = 12, height = 18)
  print(gg1)
  dev.off()
  
  # Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
  
  # Selecting number of patterns to infer. If selectK has been done and number of K is set, then set run.selectK to false
  # and add the number of K's to the vector in the beginning 
  if (run.selectK == TRUE) {
    CC.selectK = selectK(CC.x, pattern = "outgoing")
    ggsave2(filename = paste0(Output.dir, "SelectK_Kmer_plot_outgoing_", CC.name, "_All_pathways.pdf"), plot = CC.selectK)
    
    CC.selectK = selectK(CC.x, pattern = "incoming")
    ggsave2(filename = paste0(Output.dir, "SelectK_Kmer_plot_incoming_", CC.name, "_All_pathways.pdf"), plot = CC.selectK)
    
  } else if (run.selectK == FALSE && is.null(CC.selectK.out) == FALSE && is.null(CC.selectK.in) == FALSE) {
    K.out = CC.selectK.out[CC.i]
    K.in = CC.selectK.in[CC.i]
    
    # Identifying outgoing communication
    CC.x <- identifyCommunicationPatterns(CC.x, pattern = "outgoing", k = K.out)
    # Generate a river plot
    gg1 = netAnalysis_river(CC.x, pattern = "outgoing", color.use = subtype.col, font.size = 5)
    ggsave2(filename = paste0(Output.dir, "River_plot_outgoing_signaling_", CC.name, "_All_pathways.pdf"), plot = gg1, height = 25, width = 7)
    gg1 = netAnalysis_river(CC.x, pattern = "outgoing", color.use = subtype.col, font.size = 5, signaling = CC.pathway_name.prop.20$Category)
    ggsave2(filename = paste0(Output.dir, "River_plot_outgoing_signaling_", CC.name, "_Top20_pathways.pdf"), plot = gg1, height = 9)
    gg1 = netAnalysis_river(CC.x, pattern = "outgoing", color.use = subtype.col, font.size = 5, signaling = CC.pathway_name.prop.30$Category)
    ggsave2(filename = paste0(Output.dir, "River_plot_outgoing_signaling_", CC.name, "_Top30_pathways.pdf"), plot = gg1, height = 9)
    gg1 = netAnalysis_river(CC.x, pattern = "outgoing", color.use = subtype.col, font.size = 5, signaling = CC.pathway_name.prop.40$Category)
    ggsave2(filename = paste0(Output.dir, "River_plot_outgoing_signaling_", CC.name, "_Top40_pathways.pdf"), plot = gg1, height = 12, width = 26)
    
    # Identifying incoming communication
    CC.x <- identifyCommunicationPatterns(CC.x, pattern = "incoming", k = K.in)
    # Generate a river plot
    gg1 = netAnalysis_river(CC.x, pattern = "incoming", color.use = subtype.col, font.size = 5)
    ggsave2(filename = paste0(Output.dir, "River_plot_incoming_signaling_", CC.name, "_All_pathways.pdf"), plot = gg1, height = 25, width = 7)
    gg1 = netAnalysis_river(CC.x, pattern = "incoming", color.use = subtype.col, font.size = 5, signaling = CC.pathway_name.prop.20$Category)
    ggsave2(filename = paste0(Output.dir, "River_plot_incoming_signaling_", CC.name, "_Top20_pathways.pdf"), plot = gg1, height = 9)
    gg1 = netAnalysis_river(CC.x, pattern = "incoming", color.use = subtype.col, font.size = 5, signaling = CC.pathway_name.prop.30$Category)
    ggsave2(filename = paste0(Output.dir, "River_plot_incoming_signaling_", CC.name, "_Top30_pathways.pdf"), plot = gg1, height = 9)
    gg1 = netAnalysis_river(CC.x, pattern = "incoming", color.use = subtype.col, font.size = 5, signaling = CC.pathway_name.prop.40$Category)
    ggsave2(filename = paste0(Output.dir, "River_plot_incoming_signaling_", CC.name, "_Top40_pathways.pdf"), plot = gg1, height = 12)
    
  }
  
  for (CC.pathway in pathways.plot) {
    
    print(CC.pathway)
    
    # Chord diagram
    gg1 = netVisual_aggregate(CC.x, signaling = CC.pathway, layout = "circle")
    pdf(file=paste0(Output.dir, "CirclePlot_", CC.name, "_", CC.pathway, ".pdf"))
    print(gg1)
    dev.off()
    
    # Heatmap
    gg1 = netVisual_heatmap(CC.x, signaling = CC.pathway, color.heatmap = "Reds")
    pdf(file=paste0(Output.dir, "Heatmap_", CC.name, "_", CC.pathway, ".pdf"))
    print(gg1)
    dev.off()
    
    # Contribution of each ligand-receptor pair to the signalling pathway
    gg1 = netAnalysis_contribution(CC.x, signaling = CC.pathway)
    pdf(file=paste0(Output.dir, "Barplot_LR_pair_", CC.name, "_", CC.pathway, ".pdf"))
    print(gg1)
    dev.off()
    
    gg1 = netAnalysis_signalingRole_network(CC.x, signaling = CC.pathway)
    pdf(file=paste0(Output.dir, "Heatmap_CentralityScore_", CC.name, "_", CC.pathway, ".pdf"), width = 9, height = 7)
    print(gg1)
    dev.off()
    
    # Signaling role analysis on the aggregated cell-cell communication network from specific pathway
    gg1 <- netAnalysis_signalingRole_scatter(CC.x, signaling = CC.pathway)
    pdf(file=paste0(Output.dir, "Scatterplot_CentralityScore_", CC.name, "_", CC.pathway, ".pdf"), width = 9, height = 7)
    print(gg1)
    dev.off()
    
  }
  
  
}

groupSize <- as.numeric(table(CC.x@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(CC.x@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
gg1 = netVisual_circle(CC.x@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



