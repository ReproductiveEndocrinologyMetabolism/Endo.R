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
celltype = "Endothelial" # Epithelium_All, Immune, Endothelial, Epithelium_CtrlvsPCOS, Stromal_uSMC
DEG.comparison = "Control_PCOS"
extra.label = ""
Project_name = paste0("Endo_All_", celltype)
GO.dir = "Data/Selected_GO/"
DEG.dir = "Data/DEG_tables/"
Norm.assay = "RNA"
Group.order = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS")
CtrlvsPCOS.flag = FALSE # TRUE when working with averaged baseline sample
heatmap.break.col = 5 # Number of controls

# If CtrlvsPCOS.flag is true, automatically adjust for this
if (CtrlvsPCOS.flag == TRUE) {
  
  # Load the GO table from data dir with baseline GO
  GO.xlsx = read.xlsx(xlsxFile = paste0(GO.dir, celltype, "_CtrlvsPCOS_Curated_GO_plotting.xlsx"))
  DEG.xlsx = read.xlsx(xlsxFile = paste0(DEG.dir, "DEG_table_all_celltypes_", DEG.comparison, ".xlsx"))
  
} else if (CtrlvsPCOS.flag == FALSE) {
  
  if (extra.label == "") {
    # Load the GO and DEG table from data dir with all groups
    GO.xlsx = read.xlsx(xlsxFile = paste0(GO.dir, celltype, "_Curated_GO_plotting.xlsx"))
    DEG.xlsx = read.xlsx(xlsxFile = paste0(DEG.dir, "DEG_table_all_celltypes_", DEG.comparison, ".xlsx"))
  } else if (extra.label != "") {
    # Load the GO table from data dir with all groups
    GO.xlsx = read.xlsx(xlsxFile = paste0(GO.dir, celltype, extra.label, "_Curated_GO_plotting.xlsx"))
    DEG.xlsx = read.xlsx(xlsxFile = paste0(DEG.dir, "DEG_table_all_celltypes_", DEG.comparison, ".xlsx"))
  }
  
}

if (dir.exists(path = paste0("Output/8_DEG_GO_plotting_", celltype)) == FALSE) {
  print(paste0("Generating output directory Output/8_DEG_GO_plotting_", celltype))
  dir.create(path = paste0("Output/8_DEG_GO_plotting_", celltype), recursive = TRUE)
  Output.dir = paste0("Output/8_DEG_GO_plotting_", celltype, "/")
} else if (dir.exists(path = paste0("Output/8_DEG_GO_plotting_", celltype)) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/8_DEG_GO_plotting_", celltype, "/")
} else {
  print("Error with output directory")
}

# Function to add suffixes to duplicated row names
add_suffix_to_duplicates = function(df) {
  # Get duplicated row names
  duplicated_names = duplicated(rownames(df))
  
  # If there are duplicates, add suffixes
  if (any(duplicated_names)) {
    dup_count = table(rownames(df))
    suffixes = ave(rownames(df), rownames(df), FUN = function(x) {
      if (dup_count[x] > 1) {
        seq_len(length(x))
      } else {
        ""
      }
    })
    #rownames(df) = paste0(rownames(df), "_", suffixes)
    rownames(df) = paste0(suffixes, "_", rownames(df))
  }
  
  return(df)
}

# Extract the averaged seurat objects
celltypes.seurat = list.dirs(Output.dir, full.names = FALSE)
celltypes.seurat = celltypes.seurat[celltypes.seurat != ""]

# Extract subtypes found in the GO table and filtered celltypes.seurat based on these. 
celltypes.GO = unique(GO.xlsx$Celltype)
celltypes.seurat = celltypes.seurat[celltypes.seurat %in% celltypes.GO]

# Extract subtypes found in the DEG.table and filtered celltypes.seurat based on these. 
celltypes.DEG = unique(GO.xlsx$Celltype)
celltypes.seurat = celltypes.seurat[celltypes.seurat %in% celltypes.GO]

# Loop over the celltypes in the GO tables
for (celltype.x in celltypes.GO) {
  
  print(paste("Plotting", celltype.x))
  
  if (CtrlvsPCOS.flag == TRUE) {
    
    # Load the average seurat object of the celltype
    subset.seurat = LoadH5Seurat(file = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_CtrlvsPCOS.h5seurat"))
    Idents(object = subset.seurat) <- "orig.ident"
    DefaultAssay(subset.seurat) = "RNA"
    subset.average = AverageExpression(subset.seurat, return.seurat = TRUE, assays = "RNA")
    subset.average = subset.average@assays$RNA@scale.data
    
    # Setting the order
    Group.order = colnames(subset.average)
    
    
  } else if (CtrlvsPCOS.flag == FALSE) {
    
    # Load the average seurat object of the celltype
    subset.average = LoadH5Seurat(file = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_average.h5seurat"))
    subset.average = subset.average@assays$RNA@scale.data
    
    # Remove heatmap.break.col as CtrlvsPCOS is not used
    heatmap.break.col = 0
    
  }
  
  # Extract the target celltype from the GO tables and a vector of the pathways
  celltype.x_GO = GO.xlsx[GO.xlsx$Celltype == celltype.x,]
  celltype.x.pathways = unique(celltype.x_GO$Description)
  
  # Extract the target celltype from the DEG table
  celltype.x_DEG = DEG.xlsx[DEG.xlsx$Celltype == celltype,]
  celltype.x_DEG = DEG.xlsx[DEG.xlsx$Subtype == celltype.x,]
  
  # Make a list to store gene expression matrices in
  pathway.matrix.list = list()
  
  # A vector to keep the heatmap row break in based on pathway gene number
  # and all genes from GO's
  heatmap.vec = c()
  geneID.vec = c()
  
  # Loop over the pathways in the specific celltype
  for (pathway.x in celltype.x.pathways) {
    
    # Extract the genes of the specific pathway and remove duplicates
    pathway.x.genes = celltype.x_GO[celltype.x_GO$Description == pathway.x,]
    pathway.x.genes = pathway.x.genes$geneID
    pathway.x.genes = strsplit(pathway.x.genes, "/")
    pathway.x.genes = unique(unlist(pathway.x.genes))
    
    # From subset.average, extract gene expression matrix of pathway genes
    pathway.x.matrix = as.data.frame(subset.average[rownames(subset.average) %in% pathway.x.genes,])
    pathway.x.matrix = pathway.x.matrix[,Group.order]
    
    if (CtrlvsPCOS.flag == FALSE) {
      pathway.x.matrix = pathway.x.matrix[order(pathway.x.matrix[["Control"]], decreasing = TRUE), ]
    }

    # Add the result to the list of matrices and number of genes to the break vector
    pathway.matrix.list[[pathway.x]] = pathway.x.matrix
    heatmap.vec = c(heatmap.vec, nrow(pathway.x.matrix))
    geneID.vec = c(geneID.vec, rownames(pathway.x.matrix))
    
  }
  
  # Combine the dataframes in pathway.matrix.list. Add suffix to duplicated geneID
  pathway.matrix <- do.call(rbind, lapply(pathway.matrix.list, add_suffix_to_duplicates))
  
  
  # Recalculate the heatmap.break.row vector
  heatmap.break.row = numeric(length(heatmap.vec) - 1)
  running.sum = heatmap.vec[1]
  for (i in 1:(length(heatmap.vec) - 1)) {
    heatmap.break.row[i] <- running.sum
    running.sum <- running.sum + heatmap.vec[i + 1]
  }
  
  # Generate a heatmap with pathways genes
  plot.x = pheatmap::pheatmap(pathway.matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                              cellwidth = 10, cellheight = 5, fontsize = 5,
                              gaps_row = heatmap.break.row,
                              gaps_col = heatmap.break.col,
                              scale = "row")
  
  plot.y = pheatmap::pheatmap(pathway.matrix, cluster_rows = FALSE, cluster_cols = TRUE,
                              cellwidth = 10, cellheight = 5, fontsize = 5,
                              gaps_row = heatmap.break.row,
                              gaps_col = heatmap.break.col,
                              scale = "row")
  
  # Generate a heatmap with all the genes from the GO table
  heatmap.genes = as.data.frame(subset.average[rownames(subset.average) %in% geneID.vec,])
  heatmap.genes = heatmap.genes[,Group.order]
  
  if (CtrlvsPCOS.flag == FALSE) {
    ggsave2(plot = plot.x, filename = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_GO_pathway_genes_heatmap_Group.pdf"))
    ggsave2(plot = plot.y, filename = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_GO_pathway_genes_heatmap_Group_clustered.pdf"))
    
    #heatmap.genes = heatmap.genes[order(heatmap.genes[["Control"]], decreasing = TRUE), ]
    
    plot.x = pheatmap::pheatmap(heatmap.genes, cluster_rows = TRUE, cluster_cols = TRUE,
                                cellwidth = 10, cellheight = 5, fontsize = 5, scale = "row")
    ggsave2(plot = plot.x, filename = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_GO_pathway_ALL_genes_heatmap_Group.pdf"))
  } else if (CtrlvsPCOS.flag == TRUE) {
    ggsave2(plot = plot.x, filename = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_GO_pathway_genes_heatmap_CtrlvsPCOS.pdf"))
    ggsave2(plot = plot.y, filename = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_GO_pathway_genes_heatmap_CtrlvsPCOS_clustered.pdf"))
    
    plot.x = pheatmap::pheatmap(heatmap.genes, cluster_rows = TRUE, cluster_cols = TRUE,
                                cellwidth = 10, cellheight = 5, fontsize = 5, scale = "row")
    ggsave2(plot = plot.x, filename = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_GO_pathway_ALL_genes_heatmap_CtrlvsPCOS_clustered.pdf"))
    
    plot.x = pheatmap::pheatmap(heatmap.genes, cluster_rows = FALSE, cluster_cols = FALSE,
                                cellwidth = 10, cellheight = 5, fontsize = 5, gaps_col = heatmap.break.col, scale = "row")
    ggsave2(plot = plot.x, filename = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_GO_pathway_ALL_genes_heatmap_CtrlvsPCOS.pdf"))
    
  }
  
  # From subset.average, extract gene expression matrix of pathway genes
  DEG.x.matrix = as.data.frame(subset.average[rownames(subset.average) %in% celltype.x_DEG$gene,])
  DEG.x.matrix = DEG.x.matrix[,Group.order]
  
  plot.z = pheatmap::pheatmap(DEG.x.matrix, cluster_rows = TRUE, cluster_cols = TRUE,
                              cellwidth = 10, cellheight = 1, fontsize = 5, show_rownames = FALSE,
                              scale = "row")
  ggsave2(plot = plot.z, filename = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_DEGs_ALL_genes_heatmap_CtrlvsPCOS_clustered.pdf"))
  
   
}

#### Loop for only doing DEG plotting

# Extract the target celltype from the DEG table
DEG.xlsx = DEG.xlsx[DEG.xlsx$Celltype == celltype,]
celltypes.DEG = unique(DEG.xlsx$Subtype)

for (celltype.x in celltypes.DEG) {
  
  print(paste("Plotting", celltype.x))
  
  # Extract only the subtype
  celltype.x_DEG = DEG.xlsx[DEG.xlsx$Subtype == celltype.x,]
  
  if (nrow(celltype.x_DEG) < 3) {
    print(paste("Less than 5 DEG's in", celltype.x))
    next()
  }
  
  # Load the average seurat object of the celltype
  subset.average = LoadH5Seurat(file = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_average.h5seurat"))
  subset.average = subset.average@assays$RNA@scale.data
  
  # From subset.average, extract gene expression matrix of pathway genes
  DEG.x.matrix = as.data.frame(subset.average[rownames(subset.average) %in% celltype.x_DEG$gene,])
  DEG.x.matrix = DEG.x.matrix[,Group.order]
  
  # Generate the heatmap and save it
  plot.z = pheatmap::pheatmap(DEG.x.matrix, cluster_rows = TRUE, cluster_cols = TRUE,
                              cellwidth = 10, cellheight = 1, fontsize = 5, show_rownames = FALSE,
                              scale = "row")
  
  ggsave2(plot = plot.z, filename = paste0(Output.dir, celltype.x, "/", Project_name, "_", celltype.x, "_DEGs_ALL_genes_heatmap_CtrlvsPCOS_clustered.pdf"))
  
  
}
