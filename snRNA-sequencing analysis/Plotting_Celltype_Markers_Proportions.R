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

celltype = "Endothelial" # Epithelium, Immune, Endothelial
Project_name = paste0("Endo_All_", celltype)
Input.dir = paste0("Output/5_", celltype, "_Clustering/")
Norm.assay = "log2"
load.average = TRUE

if (dir.exists(path = paste0("Output/Paper_plots_", celltype)) == FALSE) {
  print(paste0("Generating output directory Output/Paper_plots_", celltype))
  dir.create(path = paste0("Output/Paper_plots_", celltype), recursive = TRUE)
  Output.dir = paste0("Output/Paper_plots_", celltype, "/")
} else if (dir.exists(path = paste0("Output/Paper_plots_", celltype)) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/Paper_plots_", celltype, "/")
} else {
  print("Error with output directory")
}

# Loading reclustered Seurat
print("Seurat object loading")
endo.integrated = LoadH5Seurat(file = paste0(Input.dir, Project_name, "_reclustered_labelled.h5seurat"))
# Setting the idents
Idents(object = endo.integrated) <- paste0(celltype, "_labelled")
DefaultAssay(endo.integrated) = "RNA"
Norm.assay = "RNA"

# Setting colours and order for the groups
Group.cols = c("#A0A0A0", "#D098B2", "#95BFE1", "#65D46E")
Group.order = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS")

# Plot the UMAP
if (celltype == "Epithelium") {
  # Colour order is 
  Dimplot.colors = c("#65C2A5", #Lumenal
                     "#D3020D", #SOX9+ LGR5+
                     "#2570B7", # SOX9+ LGR5-
                     "#FAA0A1", # SOX9+ proliferative
                     "#FF7F00", # AR+
                     "#DD8D61") # Ciliated #DD8D61
  
  Figure_curated = c("PTGS1", "VTCN1", "SLC26A7",
                     "LGR5", "KRT5", "WNT7A",
                     "CPM", "IHH", "EMID1", "PPARG",
                     "C2CD4A", "SLC18A2", "PAEP", "CXCL14", # "HEY1", "SCGB2A2", "SPP1", "GPX3", "DPP4",
                     "MKI67", "HMGB2","AR","CDC20B", "CCNO", 
                     "HES6") # "FOXJ1", "PIFO", "TP73"
  heatmap.order = c("Lumenal", "SOX9+ LGR5+", 
                    "SOX9+ LGR5-", "SOX9+ proliferative", 
                    "AR+", "Ciliated")
  
  # Ordering the object
  Idents(endo.integrated) <- factor(Idents(endo.integrated), levels = heatmap.order)
  
} else if (celltype == "Immune") {
  Dimplot.colors = c("darkred", "orangered", # T-cells
                     "dodgerblue", "skyblue", "navyblue", # uNKs
                     "plum", "maroon", # uMs
                     "seagreen", "limegreen", #mDC
                     "salmon", #Treg
                     "violetred", # ILC3
                     "sandybrown", "darkkhaki", # B-cells, mast cells
                     "palegoldenrod", "rosybrown", # migratory DC, pDC
                     "dimgrey") # Undefined
  
  Dimplot.colors = c("salmon", #Treg
                     "orangered", "darkred", # T-cells
                     "violetred", # ILC3
                     "sandybrown", # B-cells
                     "dodgerblue", "skyblue", "navyblue", # uNKs
                     "palegoldenrod", "rosybrown", # pDC and migratory DC
                     "seagreen", "limegreen", # DC
                     "plum", "maroon", # uMs
                      "darkkhaki", #  mast cells
                     "dimgrey") # Undefined
  
  Figure_curated = c("FOXP3", "CD3G", "CD27", "IKZF4", "CD8A", "CCL5", "CD4", 
                                "IL7R", "RORC", #Treg, T-cell CD8, T-cell CD4, ILC3
                                "MS4A1", "IGHM", #B-celler
                                "NCAM1", "ITGA1", "SPINK2", "CSF1", "CD160", "GNLY", # uNK 1-3
                                "IL3RA", "LILRA4", "PLD4",# pDC
                                "EBI3", "CCR7", "CCL19", # Migratory DC
                                "BATF3", "CADM1", "CLEC9A", #DC1
                                "CLEC10A", "FCER1A", "CD1C", #DC2
                                "CD14", "SELENOP", "HMOX1", "IL1B", #uM 1-2
                                "CPA3", "KIT", "MS4A2") # Mast cells
  
  heatmap.order = c("Tregs", "T-cells CD8+", 
                             "T-cells CD4+", "ILC3", "B-cells",
                             "uNK 1", "uNK 2", "uNK 3",
                             "pDC", "Migratory DC", 
                             "DC1", "DC2", "uM 1", "uM 2",
                             "Mast cells", "Undefined")
  
  # Ordering the object
  Idents(endo.integrated) <- factor(Idents(endo.integrated), levels = heatmap.order)
  
  # Remove undefined cluster from immune data
  #Clean.subset = subset(endo.integrated, idents = c("Stromal", "Epithelium",
  #                                                  "Immune", "uSMC",
  #                                                  "Endothelial", "Lymphatic"))
  
} else if (celltype == "Endothelial") {
  
  # Colour order is 
  Dimplot.colors = c("#8491B4CC", # Endothelial Vein
                     "#B09C8599", # Endothelial Artery
                     "#F39B7F99", # Endothelial proliferative
                     "#00A08799", # Mesenchymal
                     "#7E6148CC") # Lymphatic
  
  Figure_curated = c("PECAM1", "CD34" ,
                     "ACKR1", "PLVAP", # Endothelial Vein
                     "CD34", "SEMA3G", "GJA5", # Endothelial Artery
                     "TOP2A", "MKI67", # Endothelial proliferative
                     "COL3A1", "WNT5A", "MMP11", # Mesenchymal
                     "PROX1", "FLT4") # Lymphatic
                     
  heatmap.order = c("Endothelial Vein", 
                    "Endothelial Artery", 
                    "Endothelial proliferative", 
                    "Mesenchymal", 
                    "Lymphatic")
  
  # Ordering the object
  Idents(endo.integrated) <- factor(Idents(endo.integrated), levels = heatmap.order)
  
} else if (celltype == "Stromal") {
  
  # Colour order is 
  Dimplot.colors = c("#E64B35CC", # Stroma 1
                     "#4DBBD599", # Stroma 2
                     "#FAA0A1", # Stroma proliferative
                     "#91D1C2CC", # Fibroblast
                     "#F39B7FCC") # uSMC
  
  Figure_curated = c("ESR1", "PGR", "IGF1", "ECM1", "PAEP", #Stroma
                     "OGN", "TOP2A", "MKI67", # Stroma proliferative
                     "THY1", "COL1A1", "PCOLCE", "C7", # Fibroblast
                     "ACTA2", "ACTG2", "MCAM") # uSMC) # Lymphatic
  
  heatmap.order = c("Stroma 1", 
                    "Stroma 2", 
                    "Stroma proliferative", 
                    "Fibroblast", 
                    "uSMC")
  
  # Ordering the object
  Idents(endo.integrated) <- factor(Idents(endo.integrated), levels = heatmap.order)
  
}

# Generating UMAPs with and wihtout labels
DimPlot(endo.integrated, cols = Dimplot.colors, raster = FALSE) + ggtitle(paste0(celltype, " UMAP, n = ", ncol(endo.integrated), " nuclei")) + NoLegend()
ggsave2(paste0(Output.dir, Project_name,"_no-labels_UMAP.pdf"), dpi = 700)

DimPlot(endo.integrated, cols = Dimplot.colors, label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir, Project_name, "_labels_UMAP_.pdf"), dpi = 700)


# Do heatmap by calculating average expression
endo.average = AverageExpression(endo.integrated, return.seurat = TRUE, assays = "RNA")
Idents(endo.average) <- factor(Idents(endo.average), levels = heatmap.order)

# Save the labelled and re-ordered object
SaveH5Seurat(endo.average, paste0(Output.dir,Project_name, "_reclustered_labelled_average.h5seurat"), overwrite = TRUE)

# Generate the heatmap
DoHeatmap(endo.average, features = Figure_curated, raster = FALSE,
          group.bar = TRUE, group.colors = c(rep("white", length(levels(endo.average)))),
          draw.lines = FALSE, angle = 45, hjust = 0) + 
  scale_fill_gradient2(low = "#2570B7", mid = "seashell", midpoint = 0, high = "#DC0000FF") +
  theme(axis.text.y = element_text(face = "italic"))

ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Marker_Heatmap_selected_genes.pdf"), 
        dpi = 700, height = 14, width = 8)

Idents(endo.average) <- factor(Idents(endo.average), levels = rev(heatmap.order))
DoHeatmap(endo.average, features = Figure_curated, raster = FALSE, 
          group.bar = TRUE, group.colors = c(rep("white", length(levels(endo.average)))),
          draw.lines = FALSE, angle = 270, hjust = 1) + 
  scale_fill_gradient2(low = "#2570B7", mid = "seashell", midpoint = 0, high = "#DC0000FF") + 
  theme(axis.text.y = element_text(face = "italic", angle = 315))

ggsave2(paste0(Output.dir, Project_name,"_", Norm.assay, "_Marker_Heatmap_flipped_selected_genes.pdf"), 
        dpi = 700, height = 14, width = 8)

# Generate a dotplot of the markers
Idents(endo.integrated) <- factor(Idents(endo.integrated), levels = heatmap.order)
DotPlot(object = endo.integrated, features = Figure_curated,
        cols = c("navajowhite", "firebrick"), dot.min = 0, dot.scale = 10, assay = "RNA") + 
  theme(axis.text.x=element_text(angle=45, face = "italic"))
ggsave2(paste0(Output.dir, Project_name, "_", Norm.assay, "_Main_cells_Markers_Dotplot_RESCALED.pdf"), 
        dpi = 700, width = 14)

# Generating tables
table_celltype = table(Idents(endo.integrated))
table_groups = table(endo.integrated$Group_Stage)
table_idents = table(endo.integrated$orig.ident)
count.table_groups = table(Idents(endo.integrated), endo.integrated$Group_Stage)
prop.table_groups = prop.table(table(Idents(endo.integrated), endo.integrated$Group_Stage), margin = 2)
count.table_idents = table(Idents(endo.integrated), endo.integrated$orig.ident)
prop.table_idents = prop.table(table(Idents(endo.integrated), endo.integrated$orig.ident), margin = 2)

# Saving the tables
saveRDS(object = table_celltype, file = paste0(Output.dir, Project_name,"_", "main_table_celltype.rds"))
saveRDS(object = table_groups, file = paste0(Output.dir, Project_name,"_", "main_table_groups.rds"))
saveRDS(object = table_idents, file = paste0(Output.dir, Project_name,"_", "main_table_idents.rds"))
saveRDS(object = count.table_groups, file = paste0(Output.dir, Project_name,"_", "count.table_groups.rds"))
saveRDS(object = prop.table_groups, file = paste0(Output.dir, Project_name,"_", "prop.table_groups.rds"))
saveRDS(object = count.table_idents, file = paste0(Output.dir, Project_name,"_", "count.table_idents.rds"))
saveRDS(object = prop.table_idents, file = paste0(Output.dir, Project_name,"_", "prop.table_idents.rds"))

# Convert to dataframe for t.test and ggplot
prop.df_groups = as.data.frame(prop.table_groups)
prop.df_idents = as.data.frame(prop.table_idents)
prop.df_idents = prop.df_idents %>%
  mutate(sample.group = case_when(
    str_detect(Var2, "^LS.*W16$") ~ "PCOS:Lifestyle",
    str_detect(Var2, "^Met.*W16$") ~ "PCOS:Metformin",
    str_detect(Var2, ".*_W0$") ~ "PCOS:Baseline",
    str_detect(Var2, "^Ctrl_.*") ~ "Control:Baseline",
  ))

# Generate barplot of proportions
ggplot(prop.df_groups, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = Dimplot.colors) +
  scale_x_discrete(limits = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                   labels = c("Control" = "Control:Baseline", "PCOS_W0" = "PCOS:Baseline", 
                              "PCOS_W16_Met" = "PCOS:Metformin", "PCOS_W16_LS" = "PCOS:Lifestyle")) +
  xlab(label = "Group") + ylab(label = "Frequency") + labs(fill = "Cell type") +
  theme_cowplot()

ggsave2(paste0(Output.dir, Project_name, "_barplot_group_celltypes_main.pdf"), dpi = 700)

ggplot(prop.df_idents, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = Dimplot.colors) + coord_flip() +
  xlab(label = "Sample") + ylab(label = "Frequency") + labs(fill = "Cell type") +
  theme_cowplot()

ggsave2(paste0(Output.dir, Project_name, "_barplot_paired_idents_celltypes_main.pdf"), dpi = 700)

idents.order = c("Ctrl_KS204", "Ctrl_KS208", "Ctrl_KS209", "Ctrl_KS210", "Ctrl_KS211",
                 "LS_KS003_W0", "LS_KS004_W0", "LS_KS027_W0", "LS_KS063_W0", 
                 "Met_KS005_W0", "Met_KS031_W0", "Met_KS047_W0", "Met_KS054_W0", "Met_KS061_W0", 
                 "Met_KS065_W0", "Met_KS068_W0", "Met_KS073_W0", 
                 "Met_KS005_W16", "Met_KS031_W16", "Met_KS047_W16", "Met_KS054_W16",
                 "Met_KS061_W16", "Met_KS065_W16", "Met_KS073_W16",
                 "LS_KS003_W16", "LS_KS004_W16", "LS_KS027_W16")

ggplot(prop.df_idents, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") + coord_flip() +
  scale_fill_manual(values = Dimplot.colors) +
  scale_x_discrete(limits = idents.order) +
  xlab(label = "Sample") + ylab(label = "Frequency") + labs(fill = "Cell type") +
  theme_cowplot()

ggsave2(paste0(Output.dir, Project_name, "_barplot_grouped_idents_celltypes_main.pdf"), dpi = 700)

# Run t.test on all celltypes of the subset
t.test_celltype <- function(celltype.x) {
  
  # Subset the specific celltype
  prop.df_idents.x = prop.df_idents[prop.df_idents$Var1 == celltype.x,]
  
  # Run Kruskal-Wallis test
  kruskal.x = kruskal.test(Freq ~ sample.group, data = prop.df_idents.x)
  
  kruskal.x$p.value
  
  if (kruskal.x$p.value < 0.05) {
    
    print(paste("Difference between groups in", celltype.x))
    
    # If kruskal shows difference, perform Pairwise Wilcoxon Rank Sum Tests between groups
    wilcox.x <- pairwise.wilcox.test(prop.df_idents.x$Freq, 
                                     prop.df_idents.x$sample.group, 
                                     p.adjust.method="BH")
    
    
    saveRDS(object = wilcox.x$p.value, file = paste0(Output.dir, Project_name,"_", celltype.x, "_wilcox_test_padj.rds"))
    
    # Testing only ctrl vs. PCOS BL
    comp.order <- list(c("Control:Baseline", "PCOS:Baseline"))
    
    plot.x = ggbarplot(prop.df_idents.x, x = "sample.group", y = "Freq",
              fill = "sample.group", palette = Group.cols, 
              order = c("Control:Baseline", "PCOS:Baseline", 
                        "PCOS:Metformin", "PCOS:Lifestyle"),
              add = "mean_sd", legend = "none", xlab = "Group", 
              ylab = "Cell type proportion", ylim = c(0,1), 
              size = 0.5, width = 0.9) +
      stat_compare_means(comparisons = comp.order, method = "wilcox.test", 
                         label = "p.format", label.y = c(0.75, 0.85, 0.95),
                         bracket.size = 0.5) + theme_cowplot() +
      theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "none")
    
    ggsave2(plot = plot.x, filename = paste0(Output.dir, Project_name,"_", celltype.x, "_barplot_wilcox.pdf"), 
            dpi = 700)
    
  } else if (kruskal.x$p.value > 0.05) {
    print(paste("No difference between groups in", celltype.x))
  }
  
}

