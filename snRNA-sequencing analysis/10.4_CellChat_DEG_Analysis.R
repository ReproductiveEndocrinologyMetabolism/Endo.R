#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(CellChat)
library(openxlsx)
library(patchwork)
library(reticulate)
library(umap)
library(ComplexHeatmap)
library(wordcloud)
library(RColorBrewer)

reticulate::use_python("/Users/gustaw.eriksson/anaconda3/bin/python", required=T)

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
FC.cut = 0.5
padj.cut = 0.05
min.cell.exp = 0.25
gsub.pattern = ".*Endo_All_CellChat_Merged_object_(.*?).rds"

order.idents = c("Lumenal", "SOX9+ LGR5+", "SOX9+ LGR5-",
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

colors.idents = c("#65C2A5", #Lumenal
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

Group.cols = c("#A0A0A0", "#D098B2", "#95BFE1", "#65D46E")

#Dotplot.cols = brewer.pal(10, "Reds")
#display.brewer.pal(10, "Reds")

# List CellChat object to load
CC.merged.list = list.files(path = Input.dir, pattern = paste0(Project_name, ".*_CellChat_Merged_object_.*.rds"), full.names = TRUE)
CC.merged.all = list.files(path = Input.dir, pattern = paste0(Project_name, "_CellChat_Merged_object_All.rds"), full.names = TRUE)
CC.merged.list = CC.merged.list[CC.merged.list != CC.merged.all]
CC.list = list.files(path = Input.dir, pattern = paste0(Project_name, "_CellChat_object_.*.rds"), full.names = TRUE)
Group.names = sapply(CC.list, function(x) {gsub(".*_CellChat_object_(.*?).rds", "\\1", x)})
Merged.names = sapply(CC.merged.list, function(x) {gsub(".*_CellChat_Merged_object_(.*?).rds", "\\1", x)})

# Load the individual CellChat object per group
CC.list = lapply(CC.list, readRDS)
names(CC.list) = Group.names

# Editing order of vectors
CC.list = CC.list[c(1, 2, 4, 3)]
Group.names = Group.names[c(1, 2, 4, 3)]
#names(Group.cols) = Group.names

# Load the merged cellchat object
CC.merged.list = lapply(CC.merged.list, readRDS)
names(CC.merged.list) = Merged.names
CC.merged.all = readRDS(CC.merged.all)

# List and load DEG tables
DEG.list = list.files(path = "Data/DEG_tables/", full.names = TRUE)
DEG.names = sapply(DEG.list, function(x) {gsub(".*DEG_table_all_celltypes_(.*?).xlsx", "\\1", x)})
DEG.list = lapply(DEG.list, read.xlsx)
DEG.list = lapply(DEG.list, as.data.frame)

# Correct elements in dataframe
replace_subtype <- function(df) {
  within(df, Subtype[Subtype == "Luminal"] <- "Lumenal")
  within(df, Subtype[Subtype == "SOX9+ LRG5+"] <- "SOX9+ LGR5+")
}
DEG.list <- lapply(DEG.list, replace_subtype)

# Add names to DEG.list
names(DEG.list) = DEG.names

# Sanity check
gg1 <- compareInteractions(CC.merged.all, show.legend = F, group = c(1:4), color.use = Group.cols)
gg2 <- compareInteractions(CC.merged.all, show.legend = F, group = c(1:4), color.use = Group.cols, measure = "weight")
gg1 + gg2

# Generate or set the output directory
if (dir.exists(path = "Output/9_CellChat/CellChat_DEG_analysis") == FALSE) {
  print("Output/9_CellChat/CellChat_DEG_analysis")
  dir.create(path = "Output/9_CellChat/CellChat_DEG_analysis", recursive = TRUE)
  Output.dir = "Output/9_CellChat/CellChat_DEG_analysis/"
} else if (dir.exists(path = "Output/9_CellChat/CellChat_DEG_analysis") == TRUE) {
  print("Directory exists")
  Output.dir = "Output/9_CellChat/CellChat_DEG_analysis/"
} else {
  print("Error with output directory")
}

# Extract predicted interactions from the datasets
CC.interactions = lapply(CC.list, subsetCommunication)

Filter_interaction_DEG <- function(interact.x, DEG.x) {
  
  # Split the receptor column
  interact.x = separate_wider_delim(interact.x, receptor, names_sep = c("_split"), delim = "_", 
                                 cols_remove = FALSE, too_many = "drop", too_few = "align_start")
  
  # Generate new column with celltype and gene merged for celltype specific filtering
  interact.x$source_ligand = paste(interact.x$source, interact.x$ligand, sep = "_")
  interact.x$target_receptor1 = paste(interact.x$target, interact.x$receptor_split1, sep = "_")
  interact.x$target_receptor2 = paste(interact.x$target, interact.x$receptor_split2, sep = "_")
  DEG.x$CT_gene = paste(DEG.x$Subtype, DEG.x$gene, sep = "_")
  
  #Filter the interaction dataframe based on DEGs
  interact.x <- interact.x %>%
    filter(source_ligand %in% DEG.x$CT_gene | target_receptor1 %in% DEG.x$CT_gene | target_receptor2 %in% DEG.x$CT_gene)
  
  # Filter the DEG list based on predicted pathways with DEG
  DEG.x <- DEG.x %>%
    filter(DEG.x$CT_gene %in% interact.x$source_ligand | DEG.x$CT_gene %in% interact.x$target_receptor1 | DEG.x$CT_gene %in% interact.x$target_receptor2)
  
  # Return a list with the filtered dataframes
  Output.x = list(Pathway = interact.x, DEG = DEG.x)
  
  return(Output.x)
  
}

# Filter the predicted pathways based on the DEGs
interact.Control = Filter_interaction_DEG(interact.x = CC.interactions$Control, DEG.x = DEG.list$Control_PCOS)
interact.PCOS = Filter_interaction_DEG(interact.x = CC.interactions$PCOS_W0, DEG.x = DEG.list$Control_PCOS)
interact.Met = Filter_interaction_DEG(interact.x = CC.interactions$PCOS_W16_Met, DEG.x = DEG.list$PCOS_Metformin)
interact.LS = Filter_interaction_DEG(interact.x = CC.interactions$PCOS_W16_LS, DEG.x = DEG.list$PCOS_Lifestyle)

# Save tables of predicted interactions with DEG
write.xlsx(interact.Control$Pathway, file = paste0(Output.dir, "Table_Control_interaction_filtered_by_DEGs.xlsx"))
write.xlsx(interact.PCOS$Pathway, file = paste0(Output.dir, "PCOS_interaction_filtered_by_DEGs.xlsx"))
write.xlsx(interact.Met$Pathway, file = paste0(Output.dir, "Metformin_interaction_filtered_by_DEGs.xlsx"))
write.xlsx(interact.LS$Pathway, file = paste0(Output.dir, "Lifestyle_interaction_filtered_by_DEGs.xlsx"))

# Save tables of DEGs in predicted pathways
write.xlsx(interact.Control$DEG, file = paste0(Output.dir, "Table_Control_DEGs_in_pathways.xlsx"))
write.xlsx(interact.PCOS$DEG, file = paste0(Output.dir, "Table_PCOS_DEGs_in_pathways.xlsx"))
write.xlsx(interact.Met$DEG, file = paste0(Output.dir, "Table_Metformin_DEGs_in_pathways.xlsx"))
write.xlsx(interact.LS$DEG, file = paste0(Output.dir, "Table_Lifestyle_DEGs_in_pathways.xlsx"))

#### EXPLORING DATA BLOCK ####
# Overlap the filtered output
DEG.filt.Control_PCOS = interact.Control$DEG[interact.Control$DEG$CT_gene %in% interact.PCOS$DEG$CT_gene,] 
DEG.filt.PCOS_Met = merge(interact.PCOS$DEG, interact.Met$DEG, by = "CT_gene")
DEG.filt.PCOS_LS = merge(interact.PCOS$DEG, interact.LS$DEG, by = "CT_gene")

# Save tables of shared DEG genes between DEG filtered predicted tables
write.xlsx(DEG.filt.Control_PCOS, file = paste0(Output.dir, "Table_Control_PCOS_SHARED_DEGs_in_DEG_filtered_pathways.xlsx"))
write.xlsx(DEG.filt.PCOS_Met, file = paste0(Output.dir, "Table_PCOS_Metformin_SHARED_DEGs_in_DEG_filtered_pathways.xlsx"))
write.xlsx(DEG.filt.PCOS_LS, file = paste0(Output.dir, "Table_PCOS_Lifestyle_SHARED_DEGs_in_DEG_filtered_pathways.xlsx"))

# Shared interaction after filtering for DEGS
interact.filt.Control_PCOS = merge(interact.Control$Pathway, interact.PCOS$Pathway, by = c("source_ligand", "target_receptor1", "target_receptor2"))
interact.filt.PCOS_Met = merge(interact.PCOS$Pathway, interact.Met$Pathway, by = c("source_ligand", "target_receptor1", "target_receptor2"))
interact.filt.PCOS_LS = merge(interact.PCOS$Pathway, interact.LS$Pathway, by = c("source_ligand", "target_receptor1", "target_receptor2"))

write.xlsx(interact.filt.Control_PCOS, file = paste0(Output.dir, "Table_Control_PCOS_SHARED_interaction_DEG_filtered.xlsx"))
write.xlsx(interact.filt.PCOS_Met, file = paste0(Output.dir, "Table_PCOS_Metformin_SHARED_interaction_DEG_filtered.xlsx"))
write.xlsx(interact.filt.PCOS_LS, file = paste0(Output.dir, "Table_PCOS_Lifestyle_SHARED_interaction_DEG_filtered.xlsx"))

# Extracting interaction after filtering that are NOT SHARED
interact.filt.Control_independent = CC.interactions$Control[!(CC.interactions$Control$pathway_name %in% CC.interactions$PCOS_W0$pathway_name), ]
interact.filt.PCOS_independent = CC.interactions$PCOS_W0[!(CC.interactions$PCOS_W0$pathway_name %in% CC.interactions$Control$pathway_name), ]

write.xlsx(interact.filt.Control_independent, file = paste0(Output.dir, "Table_Control_interaction_NOT_SHARED_WITH_PCOS.xlsx"))
write.xlsx(interact.filt.PCOS_independent, file = paste0(Output.dir, "Table_PCOS_interaction_NOT_SHARED_WITH_Control.xlsx"))

# Check the number of pathways
table(interact.Control$Pathway$pathway_name)
table(interact.PCOS$Pathway$pathway_name)
table(interact.Met$Pathway$pathway_name)
table(interact.LS$Pathway$pathway_name)

# Run netAnalysis_computeCentrality on all cellchat objects before plotting
CC.list = lapply(CC.list, netAnalysis_computeCentrality)

# Plot out pathways found after filtering for DEGs
p.union <- union(CC.list[[1]]@netP$pathways, CC.list[[2]]@netP$pathways)
p.union = p.union[p.union %in% interact.filt.Control_PCOS$pathway_name.x]
# Extract top 30
p.union = p.union[1:30]

ht1 = netAnalysis_signalingRole_heatmap(CC.list[[1]], pattern = "outgoing", signaling = p.union, 
                                        title = names(CC.list)[1], width = 8, height = 14, color.use = colors.idents)
ht2 = netAnalysis_signalingRole_heatmap(CC.list[[2]], pattern = "outgoing", signaling = p.union, 
                                        title = names(CC.list)[2], width = 8, height = 14, color.use = colors.idents)
pdf(file=paste0(Output.dir, "Outgoing_signal_heatmap_Control_PCOS_filtered_pathways.pdf"), height = 12, width = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(CC.list[[1]], pattern = "incoming", signaling = p.union, 
                                        title = names(CC.list)[1], width = 8, height = 14, color.use = colors.idents)
ht2 = netAnalysis_signalingRole_heatmap(CC.list[[2]], pattern = "incoming", signaling = p.union, 
                                        title = names(CC.list)[2], width = 8, height = 14, color.use = colors.idents)
pdf(file=paste0(Output.dir, "Incoming_signal_heatmap_Control_PCOS_filtered_pathways.pdf"), height = 12, width = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(CC.list[[1]], pattern = "all", signaling = p.union, 
                                        title = names(CC.list)[1], width = 8, height = 14, color.use = colors.idents)
ht2 = netAnalysis_signalingRole_heatmap(CC.list[[2]], pattern = "all", signaling = p.union, 
                                        title = names(CC.list)[2], width = 8, height = 14, color.use = colors.idents)
pdf(file=paste0(Output.dir, "All_signal_heatmap_Control_PCOS_filtered_pathways.pdf"), height = 12, width = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

gg1 <- rankNet(CC.merged.all, comparison = c(1,2), mode = "comparison", stacked = T, do.stat = TRUE, color.use = Group.cols[c(1,2)], 
               signaling = p.union)

# Pathway target is required for this to work
pdf(file=paste0(Output.dir, "Control_shared_signaling_network.pdf"), height = 8, width = 8)
netAnalysis_signalingRole_network(CC.list$Control, signaling = p.union, width = 13, height = 2.5, font.size = 10, color.use = colors.idents)
dev.off()

pdf(file=paste0(Output.dir, "PCOS_shared_signaling_network.pdf"), height = 8, width = 8)
netAnalysis_signalingRole_network(CC.list$PCOS_W0, signaling = pathway.target, width = 13, height = 2.5, font.size = 10, color.use = colors.idents)
dev.off()


# Function to plot signalling of selected hormones
Paired_Multiple_heatmap <- function(CC.x = "Control", CC.y = "PCOS_W0", pathway.x, save.name, width.x = 12, height.x = 8, CC.merged = CC.merged.list$Control_PCOS) {
  
  # Skip this part when working on input with no merged CC, for example with independent pathways
  if (is.null(CC.merged) == FALSE) {
    # Generate the plots
    ht.x.out = netAnalysis_signalingRole_heatmap(CC.list[[CC.x]], signaling = pathway.x, pattern = "outgoing", color.use = colors.idents)
    ht.x.in = netAnalysis_signalingRole_heatmap(CC.list[[CC.x]], signaling = pathway.x, pattern = "incoming", color.use = colors.idents)
    ht.y.out = netAnalysis_signalingRole_heatmap(CC.list[[CC.y]], signaling = pathway.x, pattern = "outgoing", color.use = colors.idents)
    ht.y.in = netAnalysis_signalingRole_heatmap(CC.list[[CC.y]], signaling = pathway.x, pattern = "incoming", color.use = colors.idents)
    ht.x = ht.x.out + ht.x.in
    ht.y = ht.y.out + ht.y.in
    
    # Save the plots
    pdf(file=paste0(Output.dir, CC.x, "_", save.name, "_signalingRole_heatmap.pdf"), height = height.x, width = width.x)
    draw(ht.x.out + ht.x.in, ht_gap = unit(0.5, "cm"))
    dev.off()
    
    pdf(file=paste0(Output.dir, CC.y, "_", save.name, "_signalingRole_heatmap.pdf"), height = height.x, width = width.x)
    draw(ht.y.out + ht.y.in, ht_gap = unit(0.5, "cm"))
    dev.off()
  }
  
  for (pathway.target in pathway.x) {
    
    print(paste("Plotting", pathway.target))
    
    
    if (pathway.target %in% CC.interactions[[CC.x]]$pathway_name && pathway.target %in% CC.interactions[[CC.y]]$pathway_name == TRUE) {
      
      # Plot signaling roles and save the plots
      pdf(file=paste0(Output.dir, CC.x, "_", save.name, "_", pathway.target, "_shared_signaling_network.pdf"), height = 8, width = 8)
      netAnalysis_signalingRole_network(CC.list[[CC.x]], signaling = pathway.target, width = 13, height = 2.5, font.size = 10, color.use = colors.idents)
      dev.off()
      pdf(file=paste0(Output.dir, CC.y, "_", save.name, "_", pathway.target, "_shared_signaling_network.pdf"), height = 8, width = 8)
      netAnalysis_signalingRole_network(CC.list[[CC.y]], signaling = pathway.target, width = 13, height = 2.5, font.size = 10, color.use = colors.idents)
      dev.off()
      
      # Extract source and targets related to the target pathway
      df.x = CC.interactions[[CC.x]]
      df.y = CC.interactions[[CC.y]]
      df.x = df.x[df.x$pathway_name == pathway.target, ]
      df.y = df.y[df.y$pathway_name == pathway.target, ]
      source.xy = as.vector(unique(c(df.x$source, df.y$source)))
      target.xy = as.vector(unique(c(df.x$target, df.y$target)))
      
      gg.x = netVisual_bubble(CC.merged, sources.use = source.xy, targets.use = target.xy, angle.x = 45, comparison = c(1,2), signaling = pathway.target,
                       remove.isolate = TRUE, title.name = paste(pathway.target, "signalling in", CC.x, "&", CC.y))
      ggsave2(plot = gg.x, filename = paste0(Output.dir, CC.x, "_", CC.y, "_", save.name, "_", pathway.target, "_shared_signaling_bubble.pdf"), height = 8, width = 8)
      
      
    } else if (pathway.target %in% CC.interactions[[CC.x]]$pathway_name && pathway.target %in% CC.interactions[[CC.y]]$pathway_name == FALSE) {
      print(paste(pathway.target, "NOT SHARED BY BOTH DATASETS"))
      
      if (pathway.target %in% CC.interactions[[CC.x]]$pathway_name) {
        # Plot signaling roles and save the plots
        print(paste(pathway.x, "in Control"))
        pdf(file=paste0(Output.dir, CC.x, "_", save.name, "_", pathway.target, "_signaling_network.pdf"), height = 8, width = 8)
        netAnalysis_signalingRole_network(CC.list[[CC.x]], signaling = pathway.target, width = 13, height = 2.5, font.size = 10, color.use = colors.idents)
        dev.off()
      } else if (pathway.target %in% CC.interactions[[CC.y]]$pathway_name) {
        # Plot signaling roles and save the plots
        pdf(file=paste0(Output.dir, CC.y, "_", save.name, "_", pathway.target, "_signaling_network.pdf"), height = 8, width = 8)
        netAnalysis_signalingRole_network(CC.list[[CC.y]], signaling = pathway.target, width = 13, height = 2.5, font.size = 10, color.use = colors.idents)
        dev.off()
      }
      
    } 
    
  }
  
}

# Plot hormone interactions, not specifically with DEGs.
CC.hormones.Control = CC.interactions$Control[CC.interactions$Control$annotation == "Non-protein Signaling",]
CC.hormones.PCOS = CC.interactions$PCOS_W0[CC.interactions$PCOS_W0$annotation == "Non-protein Signaling",]
CC.hormones.pathways = as.vector(unique(c(CC.hormones.Control$pathway_name, CC.hormones.PCOS$pathway_name)))
plot.hormone.signalling = Paired_Multiple_heatmap(pathway.x = CC.hormones.pathways, save.name = "Hormones")

# Function to extract source, target and interactions from pathway
Pathway_extractor <- function(pathway.x, CC.x) {
  
  # Get the columns name in case the dataframes as rbinded
  interaction.col = names(CC.x)[grep("interaction_name", names(CC.x))][1]
  pathway.col = names(CC.x)[grep("pathway_name", names(CC.x))][1]
  
  if (ncol(CC.x) > 11) {
    source.col = names(CC.x)[grep("source", names(CC.x))][2]
    target.col = names(CC.x)[grep("target", names(CC.x))][3]
    ligand.col = names(CC.x)[grep("ligand", names(CC.x))][2]
    receptor.col = names(CC.x)[grep("receptor", names(CC.x))][c(3,4)]
  } else if (ncol(CC.x) <= 11) {
    source.col = names(CC.x)[grep("source", names(CC.x))][1]
    target.col = names(CC.x)[grep("target", names(CC.x))][1]
    ligand.col = names(CC.x)[grep("ligand", names(CC.x))][1]
    receptor.col = names(CC.x)[grep("receptor", names(CC.x))][1]
  }
  
  
  # Filter the CellChat dataframe
  CC.x = CC.x[CC.x[pathway.col] == pathway.x,]
  
  # Remove any levels
  CC.x = droplevels(CC.x)
  
  # Extract interactions, sources and pathways, and output the data
  output.list = list()
  output.list$interaction = as.vector(unique(CC.x[[interaction.col]]))
  output.list$source = as.vector(unique(CC.x[[source.col]]))
  output.list$target = as.vector(unique(CC.x[[target.col]]))
  output.list$ligand = as.vector(unique(CC.x[[ligand.col]]))
  output.list$receptor = unique(c(as.vector(unique(CC.x[[receptor.col[1]]])), 
                                  as.vector(unique(CC.x[[receptor.col[2]]]))))
  
  # If the input dataframe has non-split ligand and receptors, then perform split
  if (any(grepl("_", output.list$ligand))) {
    output.list$ligand = unique(unlist(strsplit(output.list$ligand, split = "_")))
  }
  if (any(grepl("_", output.list$receptor))) {
    output.list$receptor = unique(unlist(strsplit(output.list$receptor, split = "_")))
  }
  
  # Make interaction and pathway as dataframe for plotting
  output.list$pair.df = as.data.frame(CC.x[interaction.col])
  colnames(output.list$pair.df) = "interaction_name"
  
  output.list = lapply(output.list, na.omit)

  return(output.list)

}

# Interesting pathways
table(interact.filt.Control_PCOS$pathway_name.x)

# Cell-cell contact
CADM.extract = Pathway_extractor(pathway.x = "CADM", CC.x = interact.filt.Control_PCOS)
SEMA6.extract = Pathway_extractor(pathway.x = "SEMA6", CC.x = interact.filt.Control_PCOS)
NRXN.extract = Pathway_extractor(pathway.x = "NRXN", CC.x = interact.filt.Control_PCOS)
CNTN.extract = Pathway_extractor(pathway.x = "CNTN", CC.x = interact.filt.Control_PCOS)
Cell_Cell.network = Paired_Multiple_heatmap(pathway.x = c("CADM", "SEMA6", "NRXN", "CNTN"), save.name = "Cell_Cell_Contact_Network")

# ECM receptors
COLLAGEN.extract = Pathway_extractor(pathway.x = "COLLAGEN", CC.x = interact.filt.Control_PCOS)
LAMININ.extract = Pathway_extractor(pathway.x = "LAMININ", CC.x = interact.filt.Control_PCOS)
FN1.extract = Pathway_extractor(pathway.x = "FN1", CC.x = interact.filt.Control_PCOS)
ECM.network = Paired_Multiple_heatmap(pathway.x = c("COLLAGEN", "LAMININ", "FN1"), save.name = "ECM_Network")

# Secreted signaling
IGF.extract = Pathway_extractor(pathway.x = "IGF", CC.x = interact.filt.Control_PCOS)
BMP.extract = Pathway_extractor(pathway.x = "BMP", CC.x = interact.filt.Control_PCOS)
SLIT.extract = Pathway_extractor(pathway.x = "SLIT", CC.x = interact.filt.Control_PCOS)
SPP1.extract = Pathway_extractor(pathway.x = "SPP1", CC.x = interact.filt.Control_PCOS)
#MK.extract = Pathway_extractor(pathway.x = "MK", CC.x = interact.filt.Control_PCOS)
SEMA3.extract = Pathway_extractor(pathway.x = "SEMA3", CC.x = interact.filt.Control_PCOS)
Secreted.network = Paired_Multiple_heatmap(pathway.x = c("IGF", "BMP", "SLIT", "SPP1", "SEMA3"), 
                                           save.name = "Secreted_signaling_Network")
# Hormone
DHEA.extract = Pathway_extractor(pathway.x = "DHEA", CC.x = interact.filt.Control_PCOS) # Hormone

# Control or PCOS specific pathways
table(interact.filt.Control_independent$pathway_name)
VCAM_Control.extract = Pathway_extractor(pathway.x = "VCAM", CC.x = interact.filt.Control_independent)
TRAIL_Control.extract = Pathway_extractor(pathway.x = "TRAIL", CC.x = interact.filt.Control_independent)
HGF_Control.extract = Pathway_extractor(pathway.x = "HGF", CC.x = interact.filt.Control_independent)
Testosterone_Control.extract = Pathway_extractor(pathway.x = "Testosterone", CC.x = interact.filt.Control_independent) # Hormone
#Control_specific.network = Paired_Multiple_heatmap(pathway.x = c("VCAM", "TRAIL", "HGF", "Testosterone"), 
#                                                   save.name = "Control_Specific_Network")

table(interact.filt.PCOS_independent$pathway_name)
DHT_PCOS.extract = Pathway_extractor(pathway.x = "DHT", CC.x = interact.filt.PCOS_independent) # Hormones
COMPLEMENT_PCOS.extract = Pathway_extractor(pathway.x = "COMPLEMENT", CC.x = interact.filt.PCOS_independent)

# Do a bubbleplot of interaction filtered by DEGs
Do.Bubbleplot = function(CC.merged = CC.merged.list$Control_PCOS, source.x, target.x, pathway.x, pair.x, save.name) {
  
  gg.x = netVisual_bubble(CC.merged, sources.use = source.x, targets.use = target.x, pairLR.use = pair.x,
                          angle.x = 45, comparison = c(1,2), remove.isolate = TRUE)
  ggsave2(plot = gg.x, filename = paste0(Output.dir, save.name, "_signaling_bubble.pdf"), height = 8, width = 8)
  
}

# Function to generate a signalling heatmap of selected pathways
Do.Signalling_heatmap <- function(pathway.x, save.name, width.x = 10, height.x = 12) {
  
  # Generate plots of outgoing signalling 
  ht1 = netAnalysis_signalingRole_heatmap(CC.list[[1]], pattern = "outgoing", signaling = pathway.x, 
                                          title = names(CC.list)[1], width = width.x, height = height.x, color.use = colors.idents)
  ht2 = netAnalysis_signalingRole_heatmap(CC.list[[2]], pattern = "outgoing", signaling =pathway.x, 
                                          title = names(CC.list)[2], width = width.x, height = height.x, color.use = colors.idents)
  pdf(file=paste0(Output.dir, "Outgoing_signal_heatmap_Control_PCOS_", save.name, ".pdf"), height = 12, width = 10)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  # Generate plots of incoming signalling 
  ht1 = netAnalysis_signalingRole_heatmap(CC.list[[1]], pattern = "incoming", signaling = pathway.x, 
                                          title = names(CC.list)[1], width = width.x, height = height.x, color.use = colors.idents)
  ht2 = netAnalysis_signalingRole_heatmap(CC.list[[2]], pattern = "incoming", signaling = pathway.x, 
                                          title = names(CC.list)[2], width = width.x, height = height.x, color.use = colors.idents)
  pdf(file=paste0(Output.dir, "Incoming_signal_heatmap_Control_PCOS_", save.name, ".pdf"), height = 12, width = 10)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
}

# Do signalling heatmaps of selected pathways
ECM.signalling_heatmap = Do.Signalling_heatmap(pathway.x = c("COLLAGEN", "LAMININ", "FN1"), save.name = "ECM_Signalling",
                                               width.x = 8, height.x = 1)
Cell_Cell.signalling_heatmap = Do.Signalling_heatmap(pathway.x = c("CADM", "CNTN", "NRXN"), save.name = "Cell_Cell_Signalling",
                                               width.x = 8, height.x = 1.3)
Secreted.signalling_heatmap = Do.Signalling_heatmap(pathway.x = c("IGF", "BMP", "SLIT", "SPP1", "SEMA3"), save.name = "Secreted_Signalling",
                                               width.x = 8, height.x = 1.6)
Hormone.signalling_heatmap = Do.Signalling_heatmap(pathway.x = "DHEA", save.name = "Hormone_Signalling",
                                                    width.x = 8, height.x = 0.33)
All_targeted.signalling_heatmap = Do.Signalling_heatmap(pathway.x = c("COLLAGEN", "LAMININ", "FN1",
                                                                      "CADM", "CNTN", "NRXN",
                                                                      "IGF", "BMP", "SLIT", "SPP1", "SEMA3",
                                                                      "DHEA"), save.name = "All_Targeted_Signalling", 
                                                        width.x = 8, height.x = 4)

# Do violinplots and dotplots of selected genes and celltypes
Do.GeneExpression.Plot <- function(CC.merged = CC.merged.list$Control_PCOS, CC.x = CC.list$Control, CC.y = CC.list$PCOS_W0,
                                   source.x, target.x, ligand.x, receptor.x, pathway.x = NULL, 
                                   save.name, do.vln = FALSE, set.dot.min = 0.25, set.height = NULL, set.width = NULL) {
  
  # Do the violinplots or dotplots
  if (do.vln == TRUE) {
    gg.ligand = plotGeneExpression(CC.merged, features = ligand.x, split.by = "datasets", idents = source.x, color.use = Group.cols[c(1,2)])
    gg.receptor = plotGeneExpression(CC.merged, features = receptor.x, split.by = "datasets", idents = source.x, color.use = Group.cols[c(1,2)])
    
    ggsave2(plot = gg.ligand, filename = paste0(Output.dir, save.name, "_Ligand_Vln_plot.pdf"), height = 6, width = 7)
    ggsave2(plot = gg.receptor, filename = paste0(Output.dir, save.name, "_Receptor_Vln_plot.pdf"), height = 6, width = 7)
    
  } else if (do.vln == FALSE) {
    
    # Plot merged Dotplots of both control and PCOS
    gg.ligand = plotGeneExpression(CC.merged, features = rev(ligand.x), split.by = "datasets", idents = source.x, 
                                   type = "dot", color.use = "RdYlBu", scale = TRUE, col.min = -1, col.max = 1, dot.min = set.dot.min,
                                   scale.by = "size") + scale_x_discrete(position = "top") + 
                                  theme(panel.border = element_rect(color = "black", size = 1.5), axis.text.x = element_text(angle = 90))
    gg.receptor = plotGeneExpression(CC.merged, features = rev(receptor.x), split.by = "datasets", idents = target.x,
                                   type = "dot", color.use = "RdYlBu", scale = TRUE, col.min = -1, col.max = 1, dot.min = set.dot.min,
                                   scale.by = "size") + 
      theme(panel.border = element_rect(color = "black", size = 1.5), axis.text.x = element_text(angle = 90))
    ggsave2(plot = gg.ligand, filename = paste0(Output.dir, save.name, "_Merged_Ligand_Dotplot.pdf"), height = set.height, width = set.width)
    ggsave2(plot = gg.receptor, filename = paste0(Output.dir, save.name, "_Merged_Receptor_Dotplot.pdf"), height = set.height, width = set.width)
    
    # Plot RankNet barplot of selected pathways of the plotted sources, target
    if (is.null(pathway.x) == FALSE) {
      gg.RankNet = rankNet(CC.merged, comparison = c(1,2), mode = "comparison", stacked = T, do.stat = TRUE, color.use = Group.cols[c(1,2)], 
              signaling = pathway.x)
      ggsave2(plot = gg.RankNet, filename = paste0(Output.dir, save.name, "_RankNet_barplot.pdf"))
      
    }
    
    # Plot only control
    gg.ligand = plotGeneExpression(CC.x, features = ligand.x, idents = source.x,
                                   type = "dot", color.use = "YlOrRd", scale = TRUE, col.min = 0, col.max = 1, dot.min = set.dot.min,
                                   scale.by = "size") + scale_x_discrete(position = "top") + 
      theme(panel.border = element_rect(color = "black", size = 1.5))
    gg.receptor = plotGeneExpression(CC.x, features = receptor.x, idents = target.x,
                                     type = "dot", color.use = "YlOrRd", scale = TRUE, col.min = 0, col.max = 1, dot.min = set.dot.min,
                                     scale.by = "size") + 
      theme(panel.border = element_rect(color = "black", size = 1.5))
    ggsave2(plot = gg.ligand, filename = paste0(Output.dir, save.name, "_Control_Ligand_Dotplot.pdf"), height = set.height, width = set.width/2)
    ggsave2(plot = gg.receptor, filename = paste0(Output.dir, save.name, "_Control_Receptor_Dotplot.pdf"), height = set.height, width = set.width/2)
    
    # Plot only PCOS
    gg.ligand = plotGeneExpression(CC.y, features = ligand.x, idents = source.x,
                                   type = "dot", color.use = "YlOrRd", scale = TRUE, col.min = 0, col.max = 1, dot.min = set.dot.min,
                                   scale.by = "size") + scale_x_discrete(position = "top") + 
      theme(panel.border = element_rect(color = "black", size = 1.5))
    gg.receptor = plotGeneExpression(CC.y, features = receptor.x, idents = target.x,
                                     type = "dot", color.use = "YlOrRd", scale = TRUE, col.min = 0, col.max = 1, dot.min = set.dot.min,
                                     scale.by = "size") + 
      theme(panel.border = element_rect(color = "black", size = 1.5))
    ggsave2(plot = gg.ligand, filename = paste0(Output.dir, save.name, "_PCOS_Ligand_Dotplot.pdf"), height = set.height, width = set.width/2)
    ggsave2(plot = gg.receptor, filename = paste0(Output.dir, save.name, "_PCOS_Receptor_Dotplot.pdf"), height = set.height, width = set.width/2)
    
  }
  
  
}

#######Cell-cell contact
CADM.Genexpression = Do.GeneExpression.Plot(source.x = CADM.extract$source, target.x = CADM.extract$target,
                                            ligand.x = CADM.extract$ligand, receptor.x = CADM.extract$receptor, 
                                            save.name = "Control_PCOS_CADM_Cell_contact")

CNTN.Genexpression = Do.GeneExpression.Plot(source.x = CNTN.extract$source, target.x = CNTN.extract$target,
                                            ligand.x = CNTN.extract$ligand, receptor.x = CNTN.extract$receptor, 
                                            save.name = "Control_PCOS_CNTN_Cell_contact")
NRXN.Genexpression = Do.GeneExpression.Plot(source.x = NRXN.extract$source, target.x = NRXN.extract$target,
                                            ligand.x = NRXN.extract$ligand, receptor.x = NRXN.extract$receptor, 
                                            save.name = "Control_PCOS_NRXN_Cell_contact")

# Merge info of cell-cell signaling pathways and plot together
source.cell_contact = unique(c(CADM.extract$source, CNTN.extract$source, NRXN.extract$source))
target.cell_contact = unique(c(CADM.extract$target, CNTN.extract$target, NRXN.extract$target))
ligand.cell_contact = unique(c(CADM.extract$ligand, CNTN.extract$ligand, NRXN.extract$ligand))
receptor.cell_contact = unique(c(CADM.extract$receptor, CNTN.extract$receptor, NRXN.extract$receptor))

# Filter for only celltypes with DEGs
cell_contact.filt.source = DEG.list$Control_PCOS[DEG.list$Control_PCOS$gene %in% ligand.cell_contact,]
cell_contact.filt.target = DEG.list$Control_PCOS[DEG.list$Control_PCOS$gene %in% receptor.cell_contact,]

cell_contact.filt.source = unique(cell_contact.filt.source$Subtype)
cell_contact.filt.target = unique(cell_contact.filt.target$Subtype)

#cell_contact.filt.source = cell_contact.filt.source[-9]
#cell_contact.filt.target = cell_contact.filt.target

# Plot the cell-cell contact pathways
Cell_Cell.receptor.order = c("CADM1", "NECTIN3", "NOTCH1", "NOTCH2", "NRCAM", "CLSTN1", "NLGN1")
Cell_Cell_Signaling.Genexpression = Do.GeneExpression.Plot(source.x = cell_contact.filt.source, 
                                                           target.x = c(cell_contact.filt.target, target.cell_contact[c(4, 12:13, 15, 21:22)]),
                                             ligand.x = ligand.cell_contact, 
                                             receptor.x = Cell_Cell.receptor.order,
                                             save.name = "Control_PCOS_ALL_FILTERED_Cell_contact",
                                             set.dot.min = 0.1, do.vln = FALSE, set.height = 5, set.width = 7.5)

Cell_Cell_Signaling.Genexpression = Do.GeneExpression.Plot(source.x = cell_contact.filt.source, target.x = cell_contact.filt.target,
                                                           ligand.x = sort(ligand.cell_contact), receptor.x = sort(receptor.cell_contact), 
                                                           save.name = "Control_PCOS_ALL_Cell_contact", pathway.x = c("CADM", "NRXN", "CNTN"),
                                                           set.dot.min = 0.1, do.vln = FALSE, set.height = 3.6, set.width = 5.5)

# Keep only celltypes with DEGs and L/R with DEGs
cell_contact.DEG = DEG.list$Control_PCOS[DEG.list$Control_PCOS$Subtype %in% c(cell_contact.filt.source, cell_contact.filt.target),]
cell_contact.DEG.ligand = cell_contact.DEG[cell_contact.DEG$Subtype %in% cell_contact.filt.source & cell_contact.DEG$gene %in% ligand.cell_contact, ]
cell_contact.DEG.receptor = cell_contact.DEG[cell_contact.DEG$Subtype %in% cell_contact.filt.target& cell_contact.DEG$gene %in% receptor.cell_contact, ]

DEG.cell_contact.Genexpression = Do.GeneExpression.Plot(source.x = unique(cell_contact.DEG.ligand$Subtype), target.x = unique(cell_contact.DEG.receptor$Subtype),
                                               ligand.x = sort(unique(cell_contact.DEG.ligand$gene)), receptor.x = sort(unique(cell_contact.DEG.receptor$gene)),
                                               save.name = "Control_PCOS_ONLY_DEG_cell_contact", set.dot.min = 0.1, do.vln = FALSE,
                                               set.height = 3, set.width = 6)

###### ECM receptor #####
COLLAGEN.Genexpression = Do.GeneExpression.Plot(source.x = COLLAGEN.extract$source, target.x = COLLAGEN.extract$target,
                                                ligand.x = COLLAGEN.extract$ligand, receptor.x = COLLAGEN.extract$receptor, 
                                                save.name = "Control_PCOS_COLLAGEN_ECM")
LAMININ.Genexpression = Do.GeneExpression.Plot(source.x = LAMININ.extract$source, target.x = LAMININ.extract$target,
                                               ligand.x = LAMININ.extract$ligand, receptor.x = LAMININ.extract$receptor, 
                                               save.name = "Control_PCOS_LAMININ_ECM")
FN1.Genexpression = Do.GeneExpression.Plot(source.x = FN1.extract$source, target.x = FN1.extract$target,
                                           ligand.x = FN1.extract$ligand, receptor.x = FN1.extract$receptor, 
                                           save.name = "Control_PCOS_FN1_ECM")

# Merge info of ECM pathways and plot together
source.ECM = unique(c(COLLAGEN.extract$source, LAMININ.extract$source, FN1.extract$source))
target.ECM = unique(c(COLLAGEN.extract$target, LAMININ.extract$target, FN1.extract$target))
ligand.ECM = unique(c(COLLAGEN.extract$ligand, LAMININ.extract$ligand, FN1.extract$ligand))
receptor.ECM = unique(c(COLLAGEN.extract$receptor, LAMININ.extract$receptor, FN1.extract$receptor))

# Filter for only celltypes with DEGs
ECM.filt.source = DEG.list$Control_PCOS[DEG.list$Control_PCOS$gene %in% ligand.ECM,]
ECM.filt.target = DEG.list$Control_PCOS[DEG.list$Control_PCOS$gene %in% receptor.ECM,]

ECM.filt.source = unique(ECM.filt.source$Subtype)
ECM.filt.target = unique(ECM.filt.target$Subtype)

# Plot the ECM pathways
ECM.Genexpression = Do.GeneExpression.Plot(source.x = ECM.filt.source, target.x = ECM.filt.target,
                                           ligand.x = sort(ligand.ECM), receptor.x = receptor.ECM,
                                           pathway.x = c("COLLAGEN", "LAMININ", "FN1"),
                                           save.name = "Control_PCOS_ALL_FILTERED_ECM", set.dot.min = 0.1, do.vln = FALSE,
                                           set.height = 6, set.width = 6)

# Keep only celltypes with DEGs and L/R with DEGs
ECM.filt.source = ECM.filt.source[-7]
ECM.filt.target = ECM.filt.target[-11]
ECM.DEG = DEG.list$Control_PCOS[DEG.list$Control_PCOS$Subtype %in% c(ECM.filt.source, ECM.filt.target),]
ECM.DEG.ligand = ECM.DEG[ECM.DEG$Subtype %in% ECM.filt.source & ECM.DEG$gene %in% ligand.ECM, ]
ECM.DEG.receptor = ECM.DEG[ECM.DEG$Subtype %in% ECM.filt.target& ECM.DEG$gene %in% receptor.ECM, ]

DEG.ECM.Genexpression = Do.GeneExpression.Plot(source.x = unique(ECM.DEG.ligand$Subtype), target.x = unique(ECM.DEG.receptor$Subtype),
                                           ligand.x = sort(unique(ECM.DEG.ligand$gene)), receptor.x = sort(unique(ECM.DEG.receptor$gene)),
                                           pathway.x = c("COLLAGEN", "LAMININ", "FN1"),
                                           save.name = "Control_PCOS_ONLY_DEG_ECM", set.dot.min = 0.1, do.vln = FALSE,
                                           set.height = 5, set.width = 6)

##### Secreted signaling ######
IGF.Genexpression = Do.GeneExpression.Plot(source.x = IGF.extract$source, target.x = IGF.extract$target,
                                           ligand.x = IGF.extract$ligand, receptor.x = IGF.extract$receptor, 
                                           save.name = "Control_PCOS_IGF_Secreted")
BMP.Genexpression = Do.GeneExpression.Plot(source.x = BMP.extract$source, target.x = BMP.extract$target,
                                           ligand.x = BMP.extract$ligand, receptor.x = BMP.extract$receptor, 
                                           save.name = "Control_PCOS_BMP_Secreted")

SLIT.Genexpression = Do.GeneExpression.Plot(source.x = SLIT.extract$source, target.x = SLIT.extract$target,
                                            ligand.x = SLIT.extract$ligand, receptor.x = SLIT.extract$receptor, 
                                            save.name = "Control_PCOS_SLIT_Secreted")
SPP1.Genexpression = Do.GeneExpression.Plot(source.x = SPP1.extract$source, target.x = SPP1.extract$target,
                                            ligand.x = SPP1.extract$ligand, receptor.x = SPP1.extract$receptor, 
                                            save.name = "Control_PCOS_SPP1_Secreted")
#MK.Genexpression = Do.GeneExpression.Plot(source.x = MK.extract$source, target.x = MK.extract$target,
#                                          ligand.x = MK.extract$ligand, receptor.x = MK.extract$receptor, 
#                                          save.name = "Control_PCOS_MK_Secreted")
SEMA3.Genexpression = Do.GeneExpression.Plot(source.x = SEMA3.extract$source, target.x = SEMA3.extract$target,
                                             ligand.x = SEMA3.extract$ligand, receptor.x = SEMA3.extract$receptor, 
                                             save.name = "Control_PCOS_SEMA3_Secreted")

# Merge info of Secreted pathways and plot together
source.Secreted = unique(c(IGF.extract$source, BMP.extract$source, SLIT.extract$source, SPP1.extract$source, SEMA3.extract$source))
target.Secreted = unique(c(IGF.extract$target, BMP.extract$target, SLIT.extract$target, SPP1.extract$target, SEMA3.extract$target))
ligand.Secreted = unique(c(IGF.extract$ligand, BMP.extract$ligand, SLIT.extract$ligand, SPP1.extract$ligand, SEMA3.extract$ligand))
receptor.Secreted = unique(c(IGF.extract$receptor, BMP.extract$receptor, SLIT.extract$receptor, SPP1.extract$receptor, SEMA3.extract$receptor))

# Filter for only celltypes with DEGs
Secreted.filt.source = DEG.list$Control_PCOS[DEG.list$Control_PCOS$gene %in% ligand.Secreted,]
Secreted.filt.target = DEG.list$Control_PCOS[DEG.list$Control_PCOS$gene %in% receptor.Secreted,]

Secreted.filt.source = unique(Secreted.filt.source$Subtype)
Secreted.filt.target = unique(Secreted.filt.target$Subtype)

#Secreted.filt.source = Secreted.filt.source[-10]
#Secreted.filt.target = Secreted.filt.target[-13]

# Plot the secreted pathways
Secreted.Genexpression = Do.GeneExpression.Plot(source.x = Secreted.filt.source, target.x = Secreted.filt.target,
                                           ligand.x = sort(ligand.Secreted), receptor.x = sort(receptor.Secreted), 
                                           save.name = "Control_PCOS_ALL_FILTERED_Secreted", set.dot.min = 0.1, do.vln = FALSE,
                                           set.height = 6, set.width = 7)

# Keep only celltypes with DEGs and L/R with DEGs in secrepeted pathways
Secreted.DEG = DEG.list$Control_PCOS[DEG.list$Control_PCOS$Subtype %in% c(Secreted.filt.source, Secreted.filt.target),]
Secreted.DEG.ligand = Secreted.DEG[Secreted.DEG$Subtype %in% Secreted.filt.source & Secreted.DEG$gene %in% ligand.Secreted, ]
Secreted.DEG.receptor = Secreted.DEG[Secreted.DEG$Subtype %in% Secreted.filt.target& Secreted.DEG$gene %in% receptor.Secreted, ]

DEG.ligand.custom = c("SPP1", "GDF7", "IGF1", "SEMA3E", "SLIT2")
DEG.Secreted.Genexpression = Do.GeneExpression.Plot(source.x = unique(Secreted.DEG.ligand$Subtype), target.x = unique(Secreted.DEG.receptor$Subtype),
                                           ligand.x = DEG.ligand.custom, receptor.x = sort(unique(Secreted.DEG.receptor$gene)),
                                           save.name = "Control_PCOS_ONLY_DEG_Secreted", set.dot.min = 0.1, do.vln = FALSE,
                                           set.height = 3.5, set.width = 6)

# Plot ECM and secreted pathways together as they share multiple receptors
ECM_Secreted.Genexpression = Do.GeneExpression.Plot(source.x = unique(c(ECM.filt.source, Secreted.filt.source)), 
                                                    target.x = unique(c(ECM.filt.target, Secreted.filt.target)),
                                           ligand.x = sort(unique(c(ligand.ECM, ligand.Secreted))),
                                           receptor.x = unique(c(receptor.ECM, receptor.Secreted)),
                                           pathway.x = NULL,
                                           save.name = "Control_PCOS_ALL_FILTERED_MERGED_ECM_SECRETED", set.dot.min = 0.1, do.vln = FALSE,
                                           set.height = 9, set.width = 8)

DEG.receptor.custom = c("ITGA1", "ITGA2", "ITGA9", "ITGAV", "ITGB1", "ITGB8",
                        "CD44", "ITGA6", "ITGB6",
                        "BMPR1B", "IGF1R", "PLXND1", "ROBO2")
ECM_Secreted.Genexpression = Do.GeneExpression.Plot(source.x = unique(c(ECM.DEG.ligand$Subtype, Secreted.DEG.ligand$Subtype)), 
                                                    target.x = unique(c(ECM.DEG.receptor$Subtype, Secreted.DEG.receptor$Subtype))[-13],
                                                    ligand.x = sort(unique(c(ECM.DEG.ligand$gene, Secreted.DEG.ligand$gene))),
                                                    receptor.x = DEG.receptor.custom,
                                                    pathway.x = NULL,
                                                    save.name = "Control_PCOS_ONLY_DEG_MERGED_ECM_SECRETED", set.dot.min = 0.1, do.vln = FALSE,
                                                    set.height = 8.25, set.width = 8)

# Plot the DHEA pathway, the only one in hormone signalling pathways
DHEA.Genexpression = Do.GeneExpression.Plot(source.x = DHEA.extract$source, target.x = DHEA.extract$target,
                                             ligand.x = DHEA.extract$ligand, receptor.x = DHEA.extract$receptor, 
                                             save.name = "Control_PCOS_DHEA_Secreted")


###### Independent pathways of Control and PCOS
# For VCAM, targets and endothelium (specifically artery) and source
VCAM.Genexpression = Do.GeneExpression.Plot(source.x = c("uM 1", "uM 2"), 
                                            target.x = c("Endothelial Vein", "Endothelial Artery"),
                                                ligand.x = VCAM_Control.extract$ligand, receptor.x = VCAM_Control.extract$receptor, 
                                                save.name = "Control_PCOS_VCAM_Indie", do.vln = FALSE, set.dot.min = 0, 
                                            set.height = 7, set.width = 7)

# For COMPLEMENT, C3 is DEG in LGR5+ and Ciliated. Targets är uM and DC where DEG are present
COMPLEMENT.Genexpression = Do.GeneExpression.Plot(source.x = order.idents[c(2, 6)], 
                                                  target.x = order.idents[c(17, 24)],
                                               ligand.x = COMPLEMENT_PCOS.extract$ligand, 
                                               receptor.x = COMPLEMENT_PCOS.extract$receptor, 
                                               save.name = "Control_PCOS_COMPLEMENT_Indie", 
                                               set.dot.min = 0,
                                               do.vln = FALSE)
##########

### Plot "reversed" gene expression after treatment using heatmap. Manual check was made to select "reversed" genes and celltypes
ligand.reversed = list(COL1A1 = c("uNK 1"), COL1A2 = c("uNK 2"), COL4A1 = c("AR+", "Stroma 2"), COL6A3 = c("AR+", "uSMC"),
                       LAMC2 = c("SOX9+ LGR5+"), SPP1 = c("uM 2"), IGF1 = "AR+", SEMA3E = "SOX9+ LGR5+", NRXN3 = "Lumenal",
                       SLIT3 = "SOX9+ proliferative")
receptor.reversed = list(ITGA1 = "uSMC", ITGA2 = c("AR+", "SOX9+ LGR5+", "SOX9+ proliferative", "Fibroblast"), ITGB8 = "AR+",
                         CD44 = c("AR+", "SOX9+ LGR5+", "SOX9+ proliferative", "uSMC", "uM 1"), ITGB6 = "SOX9+ LGR5+",
                         IGF1R = "Fibroblast", ROBO2 = c("AR+", "SOX9+ LGR5+", "Fibroblast", "Stroma 1", "Stroma 2", "Stroma proliferative"),
                         CADM1 = c("Lumenal", "SOX9+ LGR5-"))


# Function to generate merged DEG list with selected DEGs based on gene and celltype
selected.list = ligand.reversed
DEG.list.x = DEG.list
selected.genes.ligand = c(sort(unique(ECM.DEG.ligand$gene)), DEG.ligand.custom, ligand.cell_contact)
selected.genes.receptor = c(DEG.receptor.custom, Cell_Cell.receptor.order)
selected.genes = selected.genes.ligand
selected.genes = selected.genes.receptor
load.seurat = average.seurat
seurat.dir = "Data/Seurat_Subset_Subtypes_Averaged/"
save.name = "Ligand_DEG"
save.name = "Receptor_DEG"
