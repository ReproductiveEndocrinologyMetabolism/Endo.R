#!/usr/bin/env Rscript
library(ggplot2)
library(cowplot)
library(dplyr)
library(openxlsx)
library(enrichR)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(stringr)

# Checking and/or generating main output dir
if (dir.exists(path = "Output/7_Extented_GO_analysis") == FALSE) {
  print("Output/7_Extented_GO_analysis")
  dir.create(path = "Output/7_Extented_GO_analysis", recursive = TRUE)
  dir.create(path = "Output/7_Extented_GO_analysis/Filtered", recursive = TRUE)
  dir.create(path = "Output/7_Extented_GO_analysis/Unfiltered", recursive = TRUE)
  Output.dir = "Output/7_Extented_GO_analysis/"
} else if (dir.exists(path = "Output/7_Extented_GO_analysis") == TRUE) {
  print("Directory exists")
  Output.dir = "Output/7_Extented_GO_analysis/"
} else {
  print("Error with output directory")
}

celltype.cluster = "Endothelial" # Epithelium, Stromal, Immune, Endothelial
set.count.cutoff = 1 # Default is 3, so 3 genes are required for the GO to be captured
GSEA.dir = paste0("Output/6_DGE_analysis_", celltype.cluster, "_labelled/") # Epithelium, Stromal, Immune, Endothelial
GO_curated = list.files(path = "Data/Selected_GO/", pattern = "_Curated_GO_plotting.xlsx", full.names = TRUE)
#GSEA.dir = "Output/6_DGE_analysis_Stromal_labelled/"
#GSEA.dir = "Output/6_DGE_analysis_Immune_labelled/"
#GSEA.dir = "Output/6_DGE_analysis_Endothelial_labelled/"

# Extract celltypes with GSEA directories
celltype.list = list.dirs(GSEA.dir, full.names = FALSE)
celltype.list = celltype.list[nchar(celltype.list) != 0]
celltype.list = celltype.list[!str_detect(celltype.list,pattern="/GSEA_enrichment_analysis")]



# Function to make merged dataframe of comparecluster BP results and KEGG results
GSEA.preprocess <- function(celltype.x, GSEA.type = "BP", count.cutoff = 3, only.duplicates = FALSE,
                            comp.types = c("CtrlvsPCOS", "PCOSvsLifestyle", "PCOSvsMetformin")) {
  
  # Load the tables based on what GSEA
  print(paste("Loading xlsx tables of", celltype.x, GSEA.type))
  if (GSEA.type == "BP") {
    GSEA.df.dir = list.files(path = paste0(GSEA.dir, celltype.x, "/GSEA_enrichment_analysis"), 
                              pattern = "compareCluster.BP.all_table.xlsx", full.names = TRUE)
    GSEA.type = "compareCluster.BP.all_table.xlsx"
  } else if (GSEA.type == "KEGG") {
    GSEA.df.dir = list.files(path = paste0(GSEA.dir, celltype.x, "/GSEA_enrichment_analysis"), 
                              pattern = "DEG_enrichKEGG.all_regulated_table.xlsx", full.names = TRUE)
    GSEA.type = "enrichKEGG.all_regulated_table.xlsx"
  }
  
  # If not GSEA dir
  if (length(GSEA.df.dir) == 0) {
    print("No GSEA")
    return()
  }
  
  # Create an empty list to store the dataframes
  GSEA.df.list <- list()
  GSEA.gene_df.list <- list()
  GSEA.df_filtered.list <- list()
  GSEA.gene_filtered_df.list <- list()
  
  # Iterate over the file paths
  for (x.path in GSEA.df.dir) {
    
    # Extract which comparison it is
    x.comp <- sub(paste0(".*?(", paste(comp.types, collapse = "|"), ").*"), "\\1", x.path)
    
    # Read the dataframe from the file
    df <- read.xlsx(x.path)
    
    # Add comparison and celltype to the dataframe
    df$Comparison = x.comp
    df$Celltype = celltype.x
    
    # Order based on adjusted p-value
    # Ordering according to -log10 adjusted p-value and if GO term is identical
    df = df[order(df$p.adjust, df$ID, decreasing = FALSE),]
    
    # Extract all the genes from the dataframe
    GSEA.genes = Extract_genes(GSEA.df = df, celltype.x = celltype.x)
    GSEA.genes$Comparison = x.comp
    
    # Filter the dataframe
    # Based on gene count in GO
    df.filtered = df[df$Count >= count.cutoff,]
    
    # Check if all is filtered away
    if (nrow(df.filtered) != 0) {
      
      # Extract all the genes from the dataframe
      GSEA.genes.filtered = Extract_genes(GSEA.df = df.filtered, celltype.x = celltype.x)
      GSEA.genes.filtered$Comparison = x.comp
      
      # Append the dataframes to the lists
      GSEA.df_filtered.list[[x.comp]] <- df.filtered
      GSEA.gene_filtered_df.list[[x.comp]] = GSEA.genes.filtered
    } else if (nrow(df.filtered) == 0) {
      print("No GO's or genes after filtering")
    }
    
    # Append the dataframes to the lists
    GSEA.df.list[[x.comp]] <- df
    GSEA.gene_df.list[[x.comp]] = GSEA.genes
    
  }
  
  # Merge the dataframes to one combined
  combined.df = do.call(rbind, GSEA.df.list)
  combined.genes_df = do.call(rbind, GSEA.gene_df.list)
  combined.df_filtered = do.call(rbind, GSEA.df_filtered.list)
  combined.genes_df_filtered = do.call(rbind, GSEA.gene_filtered_df.list)
  
  # Filter to keep only duplicated GO terms on not
  comps.n = unique(combined.df_filtered$Comparison)
  if (only.duplicates == TRUE && length(comps.n) > 1) {
    GO.dups <- duplicated(combined.df_filtered$ID) | duplicated(combined.df_filtered$ID, fromLast = TRUE)
    combined.df_filtered = combined.df_filtered[GO.dups,]
  }
  
  # Save the combined dataframe
  write.xlsx(combined.df, file = paste0(Output.dir, "Unfiltered/", celltype.x, "_Combined_",  GSEA.type))
  write.xlsx(combined.genes_df, file = paste0(Output.dir, "Unfiltered/", celltype.x, "_GeneList_",  GSEA.type))
  write.xlsx(combined.df_filtered, file = paste0(Output.dir, "Filtered/", celltype.x, "_Combined_Filtered_",  GSEA.type))
  write.xlsx(combined.genes_df_filtered, file = paste0(Output.dir, "Filtered/", celltype.x, "_GeneList_Filtered_",  GSEA.type))
  
  if (nrow(combined.df_filtered) != 0) {
    return(list(Unfiltered = combined.df, Filtered = combined.df_filtered))
  } else if (nrow(combined.df_filtered) == 0) {
    return(list(Unfiltered = combined.df))
  }
  
}

Extract_genes <- function(GSEA.df, celltype.x) {
  
  # Extract all the genes and label with up or downregulated
  df.up = GSEA.df[GSEA.df$Cluster == "Upregulated",]
  df.down = GSEA.df[GSEA.df$Cluster == "Downregulated",]
  
  # Extracting the genes from the dataframe column
  if (nrow(df.up) != 0) {
    genes.up = df.up$geneID
    genes.up = unlist(strsplit(genes.up, "/"))
    genes.up = unique(genes.up)
    df.up = data.frame(Regulation = "Upregulated", GeneID = genes.up, Celltype = celltype.x)
  }
  
  if (nrow(df.down) != 0) {
    genes.down = df.down$geneID
    genes.down = unlist(strsplit(genes.down, "/"))
    genes.down = unique(genes.down)
    df.down = data.frame(Regulation = "Downregulated", GeneID = genes.down, Celltype = celltype.x)
  }
  
  combined.genes <- rbind(df.up, df.down)
  
  return(combined.genes)
  
  
}

# Function to merge all celltypes and filter for duplicates
Combine.celltypes <- function(combined.df.list, GSEA.type = "BP", 
                              celltype.cluster) {
  
  # Create an empty dataframe to store filtered and unfiltered data
  filtered.df <- data.frame()
  unfiltered.df <- data.frame()
  
  # Iterate over each element in the combined.df.list
  for (x in combined.df.list) {
    # Check if the element has "Filtered" and "Unfiltered" sub-elements
    if ("Filtered" %in% names(x) && "Unfiltered" %in% names(x)) {
      # If both sub-elements exist, append them to the respective dataframes
      filtered.df <- rbind(filtered.df, x$Filtered)
      unfiltered.df <- rbind(unfiltered.df, x$Unfiltered)
      
    }
  }
  
  # Keep only duplicate terms
  GO.dups <- duplicated(filtered.df$ID) | duplicated(filtered.df$ID, fromLast = TRUE)
  filtered_dups.df = filtered.df[GO.dups,]
  
  # Save the tables
  write.xlsx(unfiltered.df, file = paste0(Output.dir, celltype.cluster, "_Unfiltered_All_GOs_",  GSEA.type, ".xlsx"))
  write.xlsx(filtered.df, file = paste0(Output.dir, celltype.cluster, "_Filtered_All_GOs_",  GSEA.type, ".xlsx"))
  write.xlsx(filtered_dups.df, file = paste0(Output.dir, celltype.cluster, "_Filtered_Only_Duplicates_All_GOs_",  GSEA.type, ".xlsx"))
  
  # Return the tables
  return(list(Unfiltered = unfiltered.df, Filtered = filtered.df, Filtered_dups = filtered_dups.df))
  
}

# Make list of combined dataframes of all celltypes
combined.df.list = sapply(celltype.list, function(i) {
  GSEA.preprocess(celltype.x = i, GSEA.type = "BP", 
                  comp.types = c("CtrlvsPCOS", "PCOSvsLifestyle", "PCOSvsMetformin"),
                  count.cutoff = set.count.cutoff, only.duplicates = FALSE)
}, simplify = FALSE)

# Geenrate dataframe of all celltypes
combined.celltypes.list = Combine.celltypes(combined.df.list = combined.df.list,
                                         celltype.cluster = celltype.cluster)

# Plot out curated GO's
Curated_GO.plotting <- function(x.df, log10.min = NULL, log10.max = NULL, name.x = "All") {
  
  # Keep only relevant columns
  x.df <- x.df[, c("Cluster", "group", "ID", "Description", "GeneRatio", "BgRatio", "pvalue",
                   "p.adjust", "qvalue", "geneID", "Count", "enriched", "Comparison", "Celltype")]
  
  # Log10 the adjusted p-value for plotting and order by it
  x.df$log10_padjust = -log10(x.df$p.adjust)
  
  # Make log10_padjust negative for rows where Cluster is "Downregulated"
  x.df$log10_padjust[x.df$Cluster == "Downregulated"] = -abs(x.df$log10_padjust[x.df$Cluster == "Downregulated"])
 
  # Get the max/min log10 to be used as maximum in the plot
  if (is.null(c(log10.min, log10.max))) {
    log10.min = round(min(x.df$log10_padjust + -0.5)) * 2
    log10.max = round(max(x.df$log10_padjust + 0.5)) * 2
  }
  
  # Get GO celltypes
  GO.celltypes = unique(x.df$Celltype)
  
  # Make list to append celltype dataframes to later rbind
  df.list = list()
  
  for (x.type in GO.celltypes) {
    
    # Extract one celltype per plot
    type.df = x.df[x.df$Celltype == x.type,]
    
    # Order by log10 adj. p-value and GO ID
    type.df = type.df[order(type.df$ID, type.df$Comparison, type.df$log10_padjust, decreasing = TRUE),]
    
    # Make log10_padjust negative for rows where Cluster is "Downregulated"
    type.df$log10_padjust[type.df$Cluster == "Downregulated"] = -abs(type.df$log10_padjust[type.df$Cluster == "Downregulated"])
    
    # Generate the barplots
    if (log10.min < 0 & log10.max > 0) {
      plot.res = ggplot(type.df, aes(x = reorder(Description, log10_padjust), y = log10_padjust, fill = Comparison)) + 
        geom_bar(stat = "identity", position = "stack", color = "black") + ylim(log10.min, log10.max) +
        scale_fill_manual(values = c("#D098B2", "#65D46E", "#95BFE1")) + # PCOS, Lifestyle, Metformin
        geom_hline(yintercept = 0, color = "black") + coord_flip() + scale_x_discrete(position = "top") + xlab("") + 
        theme_cowplot()
    } else if (log10.min > 0) {
      plot.res = ggplot(type.df, aes(x = reorder(Description, log10_padjust), y = log10_padjust, fill = Comparison)) + 
        geom_bar(stat = "identity", position = "stack", color = "black") + ylim(-5,log10.max) +
        scale_fill_manual(values = c("#D098B2", "#65D46E", "#95BFE1")) + # PCOS, Lifestyle, Metformin
        geom_hline(yintercept = 0, color = "black") + coord_flip() + scale_x_discrete(position = "top") + xlab("") + 
        theme_cowplot()
    }
    
    ggsave2(plot = plot.res, filename = paste0(Output.dir, "enrichGO_", name.x, "_", x.type,"_log10_padjust_barplot.pdf"), 
            dpi = 700, width = 10, height = 2)
    
    # Add a celltype suffix to the description and X/Y if duplicated
    type.df$Description = paste0(type.df$Description, "_", x.type)
    
    # Order by log10 adjusted p-value
    type.df = type.df[order(type.df$log10_padjust, decreasing = FALSE),]
    
    df.list[[x.type]] = type.df
    
  }
  
  # Combine the dataframes and factorise for plotting
  combined_df <- do.call(rbind, df.list)
  
  # Find duplicate levels in the Description column
  duplicated_levels = duplicated(combined_df$Description)
  
  # Remove duplicates and set unique levels
  unique_levels = unique(combined_df$Description)
  combined_df$Description = factor(combined_df$Description, levels = unique_levels)

  
  plot.res = ggplot(combined_df, aes(x = Description, y = log10_padjust, fill = Comparison)) + 
    geom_bar(stat = "identity", position = "stack", color = "black") + ylim(log10.min, log10.max) +
    scale_fill_manual(values = c("#D098B2", "#65D46E", "#95BFE1")) + #Ctrl vs. POCS, PCOS vs. LS, PCOS vs. Met mix of colours
    geom_hline(yintercept = 0, color = "black") + coord_flip() + scale_x_discrete(position = "top") + xlab("") + 
    theme_cowplot(font_size = 5)
  
  ggsave2(plot = plot.res, filename = paste0(Output.dir, "enrichGO_", name.x, "_all_celltypes_log10_padjust_barplot.pdf"), 
          dpi = 700, width = 4, height = 3)
  
}
read.xlsx()

GO.curated = lapply(GO_curated, read.xlsx)
Epithelium.order = c("")
Epithelium.GO_plotting = Curated_GO.plotting(x.df = GO.curated[[3]], log10.min = -15, log10.max = 5, name.x = "Epithelium_All")
Epithelium.GO_plotting = Curated_GO.plotting(x.df = GO.curated[[2]], log10.min = -5, log10.max = 5, name.x = "Epithelium_CtrlvsPCOS")
Immune.GO_plotting = Curated_GO.plotting(x.df = GO.curated[[4]], log10.min = -10, log10.max = 5, name.x = "Immune_All")
Stroma.GO_plotting = Curated_GO.plotting(x.df = GO.curated[[5]], log10.min = -20, log10.max = 5, name.x = "Stroma_All")
Endothelium.GO_plotting = Curated_GO.plotting(x.df = GO.curated[[1]], log10.min = -3, log10.max = 3, name.x = "Endothelium_CtrlvsPCOS")

