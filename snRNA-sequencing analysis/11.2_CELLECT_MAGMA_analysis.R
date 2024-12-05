#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# Set variables
project.name = "CELLECT_Combined"
set.p_value.cutoff = 0.05
input.path.list = list.files(path = "Data/CELLECT_tables/", full.names = TRUE)

CT.order.unfiltered = c("Lumenal", "SOX9p_LGR5p", "SOX9p_LGR5n", "SOX9p_prolif", "AR", "Ciliated", # EPithelium
                        "Stroma_1", "Stroma_2", "Stroma_proliferative", "Fibroblast", "uSMC", # Stromal
                        "uM_1", "uM_2", "DC1", "DC2", "pDC", "Migratory_DC", "Mast_cells", # Myeloid
                        "uNK_1", "uNK_2", "uNK_3", "T-cells_CD4", "T-cells_CD8", "Tregs", "ILC3",  "B-cells",  # Lymphoid
                        "Endothelial_Artery", "Endothelial_Vein", "Endothelial_proliferative", "Lymphatic", "Mesenchymal") # Endothelium and lymphatic

# Setting the order of groups and gwas
Group.order = c("Baseline", "Control", "PCOS")
GWAS.order = c("PCOS", "ECancer", "Endo", "T2D", "BioT", "2hG", "FI", "BMI", "WHR")

# Setting the colors
Group.colors = c("Control" = "#A0A0A0",
                 "PCOS" = "#D098B2")

display.brewer.pal(n = 10, name = "Set3")
brewer.pal(n = 10, name = "Set3")

GWAS.color = c("PCOS", 
               "ECancer", 
               "Endo", 
               "T2D", 
               "BioT", 
               "2hG", 
               "FI", 
               "BMI", 
               "WHR")

CT.colors = c("Lumenal" = "#65C2A5", 
              "SOX9p_LGR5p" = "#D3020D", 
              "SOX9p_LGR5n" = "#2570B7", 
              "SOX9p_prolif" = "#FAA0A1", 
              "AR" = "#FF7F00", 
              "Ciliated" = "#DD8D61", # EPithelium
              "Stroma_1" = "#E64B35CC", 
              "Stroma_2" = "#4DBBD599", 
              "Stroma_proliferative" = "#FAA0A6", 
              "Fibroblast" = "#91D1C2CC", 
              "uSMC" = "#F39B7FCC", # Stromal
              "uM_1" = "plum", 
              "uM_2" = "maroon", 
              "DC1" = "seagreen", 
              "DC2" = "limegreen", 
              "pDC" = "palegoldenrod", 
              "Migratory_DC" = "rosybrown", 
              "Mast_cells" = "darkkhaki", # Myeloid
              "uNK_1" = "dodgerblue", 
              "uNK_2" = "skyblue", 
              "uNK_3" = "navyblue", 
              "T-cells_CD4" = "orangered", 
              "T-cells_CD8" = "darkred", 
              "Tregs" = "salmon", 
              "ILC3" = "violetred",  
              "B-cells" = "sandybrown",  # Lymphoid
              "Endothelial_Artery" = "#B09C8599", 
              "Endothelial_Vein" = "#8491B4CC", 
              "Endothelial_proliferative" = "#F39B7F99", 
              "Lymphatic" = "#7E6148CC", 
              "Mesenchymal" = "#00A08799") # Endothelium and lymphatic

# Remove unneccasary tables
#input.path.list = input.path.list[-c(1, 3:4)]
input.path.list = input.path.list[2]

# Generate or set the output directory
if (dir.exists(path = "Output/CELLECT_plotting") == FALSE) {
  print("Output/CELLECT_plotting")
  dir.create(path = "Output/CELLECT_plotting", recursive = TRUE)
  Output.dir = "Output/CELLECT_plotting/"
} else if (dir.exists(path = "Output/CELLECT_plotting") == TRUE) {
  print("Directory exists")
  Output.dir = "Output/CELLECT_plotting/"
} else {
  print("Error with output directory")
}

# Function to read csv to then merge to one table
CELLECT_load <- function(input.path, df.sep = ";") {
  
  # Load the CELLECT table
  print("Load the matrix")
  cellect.df = read.csv(file = input.path, sep = df.sep)
  
  # Run -log10 on the p-value for downstream plotting
  cellect.df$pvalue_log10 = -log10(cellect.df$pvalue)
  
  return(cellect.df)
  
}

# Function to filter the CELLECT table to output only significant enrichments 
CELLECT_filtering <- function(input.df, p.cutoff = set.p_value.cutoff) {
  
  print("Filtering the matrix")
  
  # Filter the table to only keep significant enrichments
  input.df = input.df[input.df$pvalue <= p.cutoff,]
  
  # Return the filtered table
  return(input.df)
  
}

# Plot merged table
CELLECT_plotting <- function(input.df, ct.colors.barplot, group.colors.barplot,
                             cutoff.line = -log10(set.p_value.cutoff), do.grouped = FALSE) {

  # Generate the barplots
  #input.df.baseline = input.df[input.df$specificity_id == "Baseline",]
  plot.a = ggplot(input.df, aes(x = gwas, y = pvalue_log10, fill = annotation)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), width = 0.7) + 
    scale_fill_manual(values = ct.colors.barplot) +
    labs(x = "Annotation and Specificity ID", y = "p-value (log10)", title = "P-value (log10) by Annotation and Specificity ID") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank()) +
    geom_hline(yintercept = cutoff.line, linetype = "dashed", color = "red", linewidth = 1)
  plot.a
  
  plot.a.alt = ggplot(input.df, aes(x = gwas, y = pvalue_log10, fill = annotation)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), width = 0.7) + 
    scale_fill_manual(values = ct.colors.barplot) +
    labs(x = "Annotation and Specificity ID", y = "p-value (log10)", title = "P-value (log10) by Annotation and Specificity ID") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank())
  plot.a.alt
  
  # Variables to generate background rectangles
  gwas_levels <- unique(input.df$gwas)
  gwas_positions <- seq_along(gwas_levels)
  
  plot.a.rect = ggplot(input.df, aes(x = gwas, y = pvalue_log10, fill = annotation)) +
    geom_rect(aes(xmin = gwas_positions - 0.5, xmax = gwas_positions + 0.5, ymin = -Inf, ymax = Inf),
              data = data.frame(gwas = gwas_levels, gwas_positions = gwas_positions),
              inherit.aes = FALSE, fill = rep(c("white", "grey90"), length.out = length(gwas_levels))) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), width = 0.7) + 
    scale_fill_manual(values = ct.colors.barplot) +
    labs(x = "Annotation and Specificity ID", y = "p-value (log10)", title = "P-value (log10) by Annotation and Specificity ID") +
    theme_cowplot() +
    geom_hline(yintercept = cutoff.line, linetype = "dashed", color = "red", linewidth = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank())
  plot.a.rect
  
  plot.a.alt.rect = ggplot(input.df, aes(x = gwas, y = pvalue_log10, fill = annotation)) +
    geom_rect(aes(xmin = gwas_positions - 0.5, xmax = gwas_positions + 0.5, ymin = -Inf, ymax = Inf),
              data = data.frame(gwas = gwas_levels, gwas_positions = gwas_positions),
              inherit.aes = FALSE, fill = rep(c("white", "grey90"), length.out = length(gwas_levels))) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), width = 0.7) + 
    scale_fill_manual(values = ct.colors.barplot) +
    labs(x = "Annotation and Specificity ID", y = "p-value (log10)", title = "P-value (log10) by Annotation and Specificity ID") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank())
  plot.a.alt.rect
  
  # Reverse the input.df for the flipped plot
  input.df.rev = input.df
  input.df.rev$annotation <- factor(input.df.rev$annotation, levels = rev(unique(input.df.rev$annotation)))
  input.df.rev$specificity_id <- factor(input.df.rev$specificity_id, levels = rev(unique(input.df.rev$specificity_id)))
  input.df.rev$gwas <- factor(input.df.rev$gwas, levels = rev(unique(input.df.rev$gwas)))
  input.df.rev = input.df.rev[order(input.df.rev$annotation, input.df.rev$specificity_id, input.df.rev$gwas), ]
  
  plot.a.tilt = ggplot(input.df.rev, aes(x = gwas, y = pvalue_log10, fill = annotation)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), width = 0.7) + 
    scale_fill_manual(values = ct.colors.barplot) +
    labs(x = "Annotation and Specificity ID", y = "p-value (log10)", title = "P-value (log10) by Annotation and Specificity ID") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank()) +
    geom_hline(yintercept = cutoff.line, linetype = "dashed", color = "red", linewidth = 1) + 
    coord_flip()
  plot.a.tilt
  
  plot.a.alt.tilt = ggplot(input.df.rev, aes(x = gwas, y = pvalue_log10, fill = annotation)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), width = 0.7) + 
    scale_fill_manual(values = ct.colors.barplot) +
    labs(x = "Annotation and Specificity ID", y = "p-value (log10)", title = "P-value (log10) by Annotation and Specificity ID") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank()) + 
    coord_flip()
  plot.a.alt.tilt
  
  plot.tilt = ggplot(input.df.rev, aes(x = interaction(specificity_id, annotation), y = pvalue_log10, fill = gwas)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), width = 0.7) + 
    labs(x = "Annotation and Specificity ID", y = "p-value (log10)", title = "P-value (log10) by Annotation and Specificity ID") +
    scale_fill_brewer(palette = "Set3") + theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank()) +
    geom_hline(yintercept = cutoff.line, linetype = "dashed", color = "red", linewidth = 1) + 
    coord_flip() # + scale_x_reverse()
  plot.tilt
  
  plot.tilt.alt = ggplot(input.df.rev, aes(x = interaction(specificity_id, annotation), y = pvalue_log10, fill = gwas)) +
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.1), width = 0.7) + 
    labs(x = "Annotation and Specificity ID", y = "p-value (log10)", title = "P-value (log10) by Annotation and Specificity ID") +
    scale_fill_brewer(palette = "Set3") + theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank()) +
    coord_flip() # + scale_x_reverse()
  plot.tilt.alt
  
  # Do plot of control and PCOS
  if (do.grouped == TRUE) {
    input.df.Ctrl_PCOS = input.df[!input.df$specificity_id == "Baseline",]
    Group.colors = c("Control" = "#A0A0A0",
                     "PCOS" = "#D098B2")
    
    plot.c = ggplot(input.df.Ctrl_PCOS, aes(x = interaction(specificity_id, annotation, gwas), y = pvalue_log10, fill = specificity_id)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = group.colors.barplot) +
      labs(x = "Annotation and Specificity ID", y = "p-value (log10)", title = "P-value (log10) by Annotation and Specificity ID") +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_blank()) +
      geom_hline(yintercept = cutoff.line, linetype = "dashed", color = "red", linewidth = 1)
    plot.c 
    
    # Add plots to output list
    output.list = list(Subcluster.plot.baseline = plot.a,
                       Subcluster.all.gwas = plot.b,
                       Group.subcluster.Ctrl_PCOS = plot.c)
    
  } else if (do.grouped == FALSE) {
    
    # Add plots to output list
    output.list = list(Subcluster.plot.baseline = plot.a,
                       Subcluster.plot.baseline.alt = plot.a.alt,
                       Subcluster.plot.baseline.rect = plot.a.rect,
                       Subcluster.plot.baseline.alt.rect = plot.a.alt.rect,
                       Subcluster.plot.baseline.tilt = plot.a.tilt,
                       Subcluster.plot.baseline.alt.tilt = plot.a.alt.tilt,
                       Subcluster.all.gwas.tilt = plot.tilt,
                       Subcluster.all.gwas.tilt.alt = plot.tilt.alt)
  }
  
  # Return the list
  return(output.list)
  
}

# Run CELLECT_load to create a list of dataframes to then merge
CELLECT.unfiltered = lapply(input.path.list, CELLECT_load)

# Merge the filtered dataframes to one and save
CELLECT.unfiltered = do.call(rbind, CELLECT.unfiltered)
writexl::write_xlsx(x = CELLECT.unfiltered, path = paste0(Output.dir, project.name, "_Unfiltered_table.xlsx"))

# Remove selected celltypes before filtering
CELLECT.unfiltered <- CELLECT.unfiltered[CELLECT.unfiltered$annotation != "Undefined", ]
CELLECT.unfiltered <- CELLECT.unfiltered[CELLECT.unfiltered$annotation != "nan", ]

# Reorder the data based on subcluster and group
CELLECT.unfiltered$annotation = factor(CELLECT.unfiltered$annotation, levels = CT.order.unfiltered)
CELLECT.unfiltered$specificity_id = factor(CELLECT.unfiltered$specificity_id, levels = Group.order)
CELLECT.unfiltered$gwas = factor(CELLECT.unfiltered$gwas, levels = GWAS.order)
CELLECT.unfiltered = CELLECT.unfiltered[order(CELLECT.unfiltered$annotation, CELLECT.unfiltered$specificity_id, CELLECT.unfiltered$gwas), ]

# Do filtering of the table and save it
CELLECT.filtered = CELLECT_filtering(input.df = CELLECT.unfiltered)
writexl::write_xlsx(x = CELLECT.filtered, path = paste0(Output.dir, project.name, "_Filtered_table.xlsx"))

# CHECK TO PERHAPS BOTH PLOT FILTERED AND UNFILTERED TO SEE THE DIFFERENCE IN ENRICHMENT
# Plot the results 
CELLECT_plots_unfiltered = CELLECT_plotting(input.df = CELLECT.unfiltered, ct.colors.barplot = CT.colors, group.colors.barplot = Group.colors)
CELLECT_plots_filtered = CELLECT_plotting(input.df = CELLECT.filtered, ct.colors.barplot = CT.colors, group.colors.barplot = Group.colors)

# Check the plots
CELLECT_plots_unfiltered$Subcluster.plot.baseline
CELLECT_plots_unfiltered$Subcluster.plot.baseline.alt
CELLECT_plots_unfiltered$Subcluster.plot.baseline.rect
CELLECT_plots_unfiltered$Subcluster.plot.baseline.alt.rect
CELLECT_plots_unfiltered$Subcluster.plot.baseline.tilt
CELLECT_plots_unfiltered$Subcluster.plot.baseline.alt.tilt
CELLECT_plots_unfiltered$Subcluster.all.gwas.tilt
CELLECT_plots_unfiltered$Subcluster.all.gwas.tilt.alt

CELLECT_plots_filtered$Subcluster.plot.baseline
CELLECT_plots_filtered$Subcluster.plot.baseline.alt
CELLECT_plots_filtered$Subcluster.plot.baseline.rect
CELLECT_plots_filtered$Subcluster.plot.baseline.alt.rect
CELLECT_plots_filtered$Subcluster.plot.baseline.tilt
CELLECT_plots_filtered$Subcluster.plot.baseline.alt.tilt
CELLECT_plots_filtered$Subcluster.all.gwas.tilt
CELLECT_plots_filtered$Subcluster.all.gwas.tilt.alt

# Manually save plots
ggsave2(filename = paste0(Output.dir, project.name, "_GWAS_barplot_unfiltered_Threshold.pdf"), 
        plot = CELLECT_plots_unfiltered$Subcluster.plot.baseline, height = 8, width = 16)
ggsave2(filename = paste0(Output.dir, project.name, "_GWAS_barplot_unfiltered_Threshold_Rectang.pdf"), 
        plot = CELLECT_plots_unfiltered$Subcluster.plot.baseline.rect, height = 4, width = 20)
ggsave2(filename = paste0(Output.dir, project.name, "_GWAS_barplot_filtered_NoThreshold.pdf"), 
        plot = CELLECT_plots_filtered$Subcluster.plot.baseline.alt, height = 4, width = 6)
ggsave2(filename = paste0(Output.dir, project.name, "_GWAS_barplot_filtered_NoThreshold_Rectang.pdf"), 
        plot = CELLECT_plots_filtered$Subcluster.plot.baseline.alt.rect)
ggsave2(filename = paste0(Output.dir, project.name, "_GWAS_barplot_filtered_NoThreshold_tiltet.pdf"), 
        plot = CELLECT_plots_filtered$Subcluster.plot.baseline.alt.tilt)
