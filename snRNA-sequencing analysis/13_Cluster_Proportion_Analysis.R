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

# Set the input seurat object and target idents
seurat.x = "Output/9_main_relabelled/Endo_All_Combined_labels.h5seurat"
set.ident = "Combined_labels" #"Combined_labels" # Main_Celltypes

if (set.ident == "Combined_labels") {
  
  Subcluster.correction.flag = TRUE
  
  order.celltypes = c("Lumenal", "SOX9+ LGR5+", 
                      "SOX9+ LGR5-", "SOX9+ proliferative", 
                      "AR+", "Ciliated", "Tregs", "T-cells CD8+", 
                      "T-cells CD4+", "ILC3", "B-cells",
                      "uNK 1", "uNK 2", "uNK 3",
                      "pDC", "Migratory DC", 
                      "DC1", "DC2", "uM 1", "uM 2",
                      "Mast cells", "Undefined",
                      "Endothelial Vein", 
                      "Endothelial Artery", 
                      "Endothelial proliferative", 
                      "Mesenchymal", 
                      "Lymphatic",
                      "Stroma 1", 
                      "Stroma 2", 
                      "Stroma proliferative", 
                      "Fibroblast", 
                      "uSMC")
  
  # Make vectors of cell subtypes to be used in specific functions
  celltypes.epithelium = order.celltypes[1:6]
  celltypes.immune = order.celltypes[7:22]
  celltypes.endothelium = order.celltypes[23:27]
  celltypes.stroma = order.celltypes[28:32]
  celltypes.list = list(celltypes.epithelium, celltypes.immune, celltypes.endothelium, celltypes.stroma)
  
  colors.idents = c("#65C2A5", #Lumenal
                    "#D3020D", #SOX9+ LGR5+
                    "#2570B7", # SOX9+ LGR5-
                    "#FAA0A1", # SOX9+ proliferative
                    "#FF7F00", # AR+
                    "#DD8D61", # Ciliated #DD8D61
                    "salmon", #Treg
                    "orangered", "darkred", # T-cells
                    "violetred", # ILC3
                    "sandybrown", # B-cells
                    "dodgerblue", "skyblue", "navyblue", # uNKs
                    "palegoldenrod", "rosybrown", # pDC and migratory DC
                    "seagreen", "limegreen", # DC
                    "plum", "maroon", # uMs
                    "darkkhaki", #  mast cells
                    "dimgrey", # Undefined
                    "#8491B4CC", # Endothelial Vein
                    "#B09C8599", # Endothelial Artery
                    "#F39B7F99", # Endothelial proliferative
                    "#00A08799", # Mesenchymal
                    "#7E6148CC", # Lymphatic
                    "#E64B35CC", # Stroma 1
                    "#4DBBD599", # Stroma 2
                    "#FAA0A6", # Stroma proliferative
                    "#91D1C2CC", # Fibroblast
                    "#F39B7FCC") # uSMC 
  
} else if (set.ident == "Main_Celltypes") {
  
  Subcluster.correction.flag = FALSE
  
  order.celltypes = c("Epithelium", "Stromal", "uSMC", "Myeloid", "Lymphoid", "Endothelial", "Lymphatic")
  colors.idents = c("#4DBBD5CC", # Epithelium BLUE RGB 228, 33, 28 # NEW ORDER TO FIT THE STORY
                    "#E64B35CC", # Stromal RED RGB 230 75 53
                    "#F39B7FCC", # uSMC RED-PINKISH RGB 228, 33, 28
                    "#91D1C2CC", # Myeloid RGB 228, 33, 28
                    "#00A087CC", # Lymphoid GREEN RGB 228, 33, 28
                    "#8491B4CC", # Endothelial Grey RGB 228, 33, 28
                    "#7E6148CC") # Lymphatic Brown RGB 228, 33, 28 
}

idents.order = c("1_Control", "2_Control", "3_Control", "4_Control", "5_Control",
                 "6_PCOS_Lifestyle_W0", "7_PCOS_Lifestyle_W0", "8_PCOS_Lifestyle_W0", "9_PCOS_EA_W0",
                 "10_PCOS_Metformin_W0", "11_PCOS_Metformin_W0",  "11_PCOS_Metformin_W0",
                 "12_PCOS_Metformin_W0",  "13_PCOS_Metformin_W0", "14_PCOS_Metformin_W0",
                 "15_PCOS_Metformin_W0", "16_PCOS_Metformin_W0",  "17_PCOS_Metformin_W0",
                 "18_PCOS_Lifestyle_W16", "19_PCOS_Lifestyle_W16", "20_PCOS_Lifestyle_W16",
                 "21_PCOS_Metformin_W16", "22_PCOS_Metformin_W16", "23_PCOS_Metformin_W16",
                 "24_PCOS_Metformin_W16", "25_PCOS_Metformin_W16", "26_PCOS_Metformin_W16",
                 "27_PCOS_Metformin_W16")


# Setting colours and order for the groups
Group.cols = c("#A0A0A0", "#D098B2", "#95BFE1", "#65D46E")
Group.order = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS")

# Do plotting or generate tables/excel TRUE or FALSE
Do.plotting = TRUE
Generate.tables = FALSE
Generate.excel = TRUE

# Generating the output directory
if (dir.exists(path = paste0("Output/11_Supplementary_Data_", set.ident)) == FALSE) {
  print(paste0("Generating output directory Output/11_Supplementary_Data_", set.ident))
  dir.create(path = paste0("Output/11_Supplementary_Data_", set.ident), recursive = TRUE)
  Output.dir = paste0("Output/11_Supplementary_Data_", set.ident, "/")
} else if (dir.exists(path = paste0("Output/11_Supplementary_Data_", set.ident)) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/11_Supplementary_Data_", set.ident, "/")
} else {
  print("Error with output directory")
}

# If Generate.tables is true, seurat is loaded and tables created
if (Generate.tables == TRUE) {
  # Loading reclustered Seurat
  print("Seurat object loading")
  seurat.x = LoadH5Seurat(file = seurat.x)
  # Setting the idents
  Idents(object = seurat.x) <- set.ident
  DefaultAssay(endo.integrated) = "RNA"
  Norm.assay = "RNA"
  
  # Generating tables
  table_celltype = table(Idents(seurat.x))
  table_groups = table(seurat.x$Group_Stage)
  table_idents = table(seurat.x$orig.ident)
  count.table_groups = table(Idents(seurat.x), seurat.x$Group_Stage)
  prop.table_groups = prop.table(table(Idents(seurat.x), seurat.x$Group_Stage), margin = 2)
  count.table_idents = table(Idents(seurat.x), seurat.x$orig.ident)
  prop.table_idents = prop.table(table(Idents(seurat.x), seurat.x$orig.ident), margin = 2)
  
  # Saving the tables
  saveRDS(object = table_celltype, file = paste0(Output.dir, set.ident,"_", "main_table_celltype.rds"))
  saveRDS(object = table_groups, file = paste0(Output.dir, set.ident,"_", "main_table_groups.rds"))
  saveRDS(object = table_idents, file = paste0(Output.dir, set.ident,"_", "main_table_idents.rds"))
  saveRDS(object = count.table_groups, file = paste0(Output.dir, set.ident,"_", "count.table_groups.rds"))
  saveRDS(object = prop.table_groups, file = paste0(Output.dir, set.ident,"_", "prop.table_groups.rds"))
  saveRDS(object = count.table_idents, file = paste0(Output.dir, set.ident,"_", "count.table_idents.rds"))
  saveRDS(object = prop.table_idents, file = paste0(Output.dir, set.ident,"_", "prop.table_idents.rds")) 
}

# If Generate.tables is true, the tables are loaded and excel tables are created
if (Generate.excel == TRUE) {
  
  # Read the tables
  table_celltype = as.data.frame(readRDS(file = paste0(Output.dir, set.ident,"_", "main_table_celltype.rds")))
  table_groups = as.data.frame(readRDS(file = paste0(Output.dir, set.ident,"_", "main_table_groups.rds")))
  table_idents = as.data.frame(readRDS(file = paste0(Output.dir, set.ident,"_", "main_table_idents.rds")))
  count.table_groups = as.data.frame(readRDS(file = paste0(Output.dir, set.ident,"_", "count.table_groups.rds")))
  prop.table_groups = as.data.frame(readRDS(file = paste0(Output.dir, set.ident,"_", "prop.table_groups.rds")))
  count.table_idents = as.data.frame(readRDS(file = paste0(Output.dir, set.ident,"_", "count.table_idents.rds")))
  prop.table_idents = as.data.frame(readRDS(file = paste0(Output.dir, set.ident,"_", "prop.table_idents.rds")))
  
  # Make list of the tables and name each list element
  table.list = list(table_celltype, table_groups, table_idents, count.table_groups, prop.table_groups, count.table_idents, prop.table_idents)
  names(table.list) = c(paste0(Output.dir, set.ident,"_", "main_table_celltype.xlsx"), paste0(Output.dir, set.ident,"_", "main_table_groups.xlsx"),
                        paste0(Output.dir, set.ident,"_", "main_table_idents.xlsx"), paste0(Output.dir, set.ident,"_", "count.table_groups.xlsx"),
                        paste0(Output.dir, set.ident,"_", "prop.table_groups.xlsx"), paste0(Output.dir, set.ident,"_", "count.table_idents.xlsx"),
                        paste0(Output.dir, set.ident,"_", "prop.table_idents.xlsx"))
  
  # Lapply to save the rds tables as excel files
  lapply(seq_along(table.list), function(i) {
    write.xlsx(table.list[[i]], names(table.list)[i])
  })
  
}

# Convert to dataframe for t.test and ggplot
prop.df_groups = as.data.frame(prop.table_groups)
prop.df_idents = as.data.frame(prop.table_idents)

# Add the raw count to the tables
prop.df_groups$Count = count.table_groups$Freq
prop.df_idents$Count = count.table_idents$Freq

# Add group information to sample table
prop.df_idents = prop.df_idents %>%
  mutate(sample.group = case_when(
    str_detect(Var2, "^LS.*W16$") ~ "PCOS:Lifestyle",
    str_detect(Var2, "^Met.*W16$") ~ "PCOS:Metformin",
    str_detect(Var2, ".*_W0$") ~ "PCOS:Baseline",
    str_detect(Var2, "^Ctrl_.*") ~ "Control:Baseline",
  ))

# Calculate the percentage from the frequency
prop.df_groups$Percentage = prop.df_groups$Freq*100
prop.df_idents$Percentage = prop.df_idents$Freq*100

# Convert Var1 (celltypes) to factor and specify level order
prop.df_groups$Var1 <- factor(prop.df_groups$Var1, levels=order.celltypes)
prop.df_idents$Var1 <- factor(prop.df_idents$Var1, levels=order.celltypes)

Stacked_barplots <- function(sample.df = prop.df_idents, group.df = prop.table_groups, 
                             stacked.cols = colors.idents, stacked.order = idents.order,
                             save.label = "main") {
  
  # Generate barplot of proportions
  plot.x = ggplot(group.df, aes(fill=Var1, y=Percentage, x=Var2)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values = stacked.cols) +
    scale_x_discrete(limits = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                     labels = c("Control" = "Control:Baseline", "PCOS_W0" = "PCOS:Baseline", 
                                "PCOS_W16_Met" = "PCOS:Metformin", "PCOS_W16_LS" = "PCOS:Lifestyle")) +
    xlab(label = "Group") + ylab(label = "Percentage") + labs(fill = "Cell type") +
    theme_cowplot()
  
  ggsave2(plot = plot.x, filename = paste0(Output.dir, set.ident, "_barplot_group_celltypes_", save.label, ".pdf"), dpi = 700)
  
  plot.y = ggplot(sample.df, aes(fill=Var1, y=Percentage, x=Var2)) + 
    geom_bar(position="stack", stat="identity") + 
    scale_fill_manual(values = stacked.cols) + coord_flip() +
    xlab(label = "Sample") + ylab(label = "Percentage") + labs(fill = "Cell type") +
    theme_cowplot()
  
  ggsave2(plot = plot.y, filename =  paste0(Output.dir, set.ident, "_barplot_paired_idents_celltypes_", save.label, ".pdf"), dpi = 700)
  
  plot.z = ggplot(sample.df, aes(fill=Var1, y=Freq, x=Var2)) + 
    geom_bar(position="stack", stat="identity") + coord_flip() +
    scale_fill_manual(values = stacked.cols) +
    scale_x_discrete(limits = idents.order) +
    xlab(label = "Sample") + ylab(label = "Frequency") + labs(fill = "Cell type") +
    theme_cowplot()
  
  ggsave2(plot = plot.z, filename =  paste0(Output.dir, set.ident, "_barplot_grouped_idents_celltypes", save.label, ".pdf"), dpi = 700)
  
  return(list(plot.x, plot.y, plot.x))
  
}

# Run t.test on all celltypes of the subset
t.test_celltype <- function(celltype.x, prop.df, save.label, set.y.max = NULL) {
  
  # Subset the specific celltype
  prop.df = prop.df[prop.df$Var1 == celltype.x,]
  
  # Shorten the sample ID
  prop.df$Var2 <- sub("^.*_(KS\\d+)_.*$", "\\1", prop.df$Var2)
  
  # Run Kruskal-Wallis test
  kruskal.x = kruskal.test(Freq ~ sample.group, data = prop.df)
  
  kruskal.x$p.value
  
  # Perform Pairwise Wilcoxon Rank Sum Tests between groups
  wilcox.x <- pairwise.wilcox.test(prop.df$Freq, 
                                   prop.df$sample.group, 
                                   p.adjust.method="none",
                                   exact = TRUE)$p.value
  wilcox.x.padj <- pairwise.wilcox.test(prop.df$Freq, 
                                   prop.df$sample.group, 
                                   p.adjust.method="BH",
                                   exact = TRUE)$p.value
  
  # Save the data in a dataframe to be returned
  wilcox.df = data.frame(expand.grid(dimnames(wilcox.x)), array(wilcox.x), array(wilcox.x.padj))
  wilcox.df$Celltype = celltype.x
  
  print(paste("p-value = ", wilcox.df[1,3], "in", celltype.x))
    
  # Testing only ctrl vs. PCOS BL
  comp.order <- list(c("Control:Baseline", "PCOS:Baseline"), 
                     c("Control:Baseline", "PCOS:Metformin"), c("Control:Baseline", "PCOS:Lifestyle"))
  
  # Change freq to Percentage
  prop.df$Percentage = prop.df$Freq*100
  
  # Get max proportion and round to single digit for Y-axis limit if not set
  if (is.null(set.y.max) == TRUE) {
    y.max = max(prop.df$Percentage) + 5
  } else if (!is.null(set.y.max) == TRUE) {
    y.max = set.y.max
  }
  
  
  plot.x = ggbarplot(prop.df, x = "sample.group", y = "Percentage",
                     fill = "sample.group", palette = Group.cols, 
                     order = c("Control:Baseline", "PCOS:Baseline", 
                               "PCOS:Metformin", "PCOS:Lifestyle"),
                     add = "mean_sd", legend = "none", xlab = "Group", 
                     ylab = "Cell type nuclei proportion (%)", ylim = c(0,y.max), 
                     size = 0.5, width = 0.9) +
    stat_compare_means(comparisons = comp.order, method = "wilcox.test", 
                       label = "p.format", label.y = y.max-2.5,
                       bracket.size = 0.5) + theme_cowplot() +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "none")
  
  plot.y = ggboxplot(prop.df, x = "sample.group", y = "Percentage",
                     fill = "sample.group", palette = Group.cols, 
                     order = c("Control:Baseline", "PCOS:Baseline", 
                               "PCOS:Metformin", "PCOS:Lifestyle"),
                     add = "dotplot", legend = "none", xlab = "Group", 
                     bxp.errorbar = FALSE,  
                     ylab = "Cell type nuclei proportion (%)", ylim = c(0,y.max), 
                     size = 0.5, width = 0.9) +
    stat_compare_means(comparisons = comp.order, method = "wilcox.test", 
                       label = "p.format",
                       bracket.size = 0.5) + theme_cowplot() +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "none")
  
  if (wilcox.df[1,3] < 0.05) {
    ggsave2(plot = plot.x, filename = paste0(Output.dir, set.ident,"_", celltype.x, "_", save.label, "_barplot_wilcox_significant.pdf"), 
            dpi = 700, width = 2.5)
    ggsave2(plot = plot.y, filename = paste0(Output.dir, set.ident,"_", celltype.x, "_", save.label, "_boxplot_wilcox_significant.pdf"), 
            dpi = 700, width = 2.5)
  } else if (wilcox.df[1,3] > 0.05) {
    ggsave2(plot = plot.x, filename = paste0(Output.dir, set.ident,"_", celltype.x, "_", save.label, "_barplot_wilcox_not_significant.pdf"), 
            dpi = 700, width = 2.5)
    ggsave2(plot = plot.y, filename = paste0(Output.dir, set.ident,"_", celltype.x, "_", save.label, "_boxplot_wilcox_not_significant.pdf"), 
            dpi = 700, width = 2.5)
  }
  
  #saveRDS(object = wilcox.x$p.value, file = paste0(Output.dir, set.ident,"_", celltype.x, "_", save.label, "_wilcox_test_padj.rds"))
  return(wilcox.df)
  
}

Paired_t.test <- function(celltype.x, prop.df) {
  
  print(celltype.x)
  
  # Subset the specific celltype
  prop.df = prop.df[prop.df$Var1 == celltype.x,]
  
  # Shorten the sample ID
  prop.df$Var2 <- sub("^.*_(KS\\d+)_.*$", "\\1", prop.df$Var2)
  
  # Perform paired sample t-test on baseline vs. treatment
  # Split lifestyle and metformin and filter rows with duplicated sample ID for paired samples
  paired.metformin.df = prop.df[prop.df$sample.group %in% c("PCOS:Baseline", "PCOS:Metformin"),]
  paired.metformin.df <- paired.metformin.df[duplicated(paired.metformin.df$Var2) | duplicated(paired.metformin.df$Var2, fromLast = TRUE), ]
  
  paired.lifestyle.df = prop.df[prop.df$sample.group %in% c("PCOS:Baseline", "PCOS:Lifestyle"),]
  paired.lifestyle.df <- paired.lifestyle.df[duplicated(paired.lifestyle.df$Var2) | duplicated(paired.lifestyle.df$Var2, fromLast = TRUE), ]
  
  # Perform Wilcoxon signed rank test on paired samples
  metformin.wilcox.x = wilcox.test(Freq ~ sample.group, data = paired.metformin.df, exact = TRUE, paired = TRUE)
  lifestyle.wilcox.x = wilcox.test(Freq ~ sample.group, data = paired.lifestyle.df, exact = TRUE, paired = TRUE)
  
  # Save the results in the dataframe to be returned
  metformin.wilcox.x = data.frame(Celltype = celltype.x, Treatment = "PCOS_Metformin", p.value = metformin.wilcox.x$p.value)
  lifestyle.wilcox.x = data.frame(Celltype = celltype.x, Treatment = "PCOS_Lifestyle", p.value = lifestyle.wilcox.x$p.value)
  wilcox.df = rbind(metformin.wilcox.x, lifestyle.wilcox.x)
  
  return(wilcox.df)
  
}

# Plotting stacked barplots of all celltypes
Stacked.all = Stacked_barplots(sample.df = prop.df_idents, group.df = prop.df_groups)

# If subcluster.correction is true, the freq will be normalised according to the main cell type proportion and not total
if (Subcluster.correction.flag == TRUE) {
  
  # Function to normalise against the nuclei count of the major cell type cluster
  Subset_proportion.table <- function(celltypes.subset, prop.df) {
    
    # Subset the table for selected celltypes
    prop.df = prop.df[prop.df$Var1 %in% celltypes.subset,]
    
    # Group the table and calculate the proportions
    prop.df <- prop.df %>%
      group_by(Var2) %>%
      mutate(Freq = (Count / sum(Count)))
    
    # Calculate the percentage from the frequency
    prop.df$Percentage = prop.df$Freq*100
    prop.df_idents$Percentage = prop.df_idents$Freq*100
    

    return(prop.df)
    
  }
  
  celltypes.prop_idents = lapply(celltypes.list, function (i) Subset_proportion.table(celltypes.subset = i, prop.df = prop.df_idents))
  celltypes.prop_groups = lapply(celltypes.list, function (i) Subset_proportion.table(celltypes.subset = i, prop.df = prop.df_groups))
  
  Stacked.epithelium = Stacked_barplots(sample.df = celltypes.prop_idents[[1]], group.df = celltypes.prop_groups[[1]], 
                                        stacked.cols = colors.idents[1:6], stacked.order = order.celltypes[1:6],
                                        save.label = "Epithelium")
  Stacked.immune = Stacked_barplots(sample.df = celltypes.prop_idents[[2]], group.df = celltypes.prop_groups[[2]], 
                                    stacked.cols = colors.idents[7:22], stacked.order = order.celltypes[7:22],
                                    save.label = "Immune")
  Stacked.endothelium = Stacked_barplots(sample.df = celltypes.prop_idents[[3]], group.df = celltypes.prop_groups[[3]], 
                                         stacked.cols = colors.idents[23:27], stacked.order = order.celltypes[23:27],
                                         save.label = "Endothelium")
  Stacked.stroma = Stacked_barplots(sample.df = celltypes.prop_idents[[4]], group.df = celltypes.prop_groups[[4]], 
                                    stacked.cols = colors.idents[28:32], stacked.order = order.celltypes[28:32],
                                    save.label = "Stroma")
  
  prop.df_idents = bind_rows(celltypes.prop_idents)
  prop.df_groups = bind_rows(celltypes.prop_groups)
  
}

print("Testing t.test")
wilcox.list = lapply(order.celltypes, function(i) t.test_celltype(celltype.x = i, prop.df = prop.df_idents, save.label = set.ident, set.y.max = 100))
lymphoid.plot = t.test_celltype(celltype.x = "Lymphoid", prop.df = prop.df_idents, save.label = paste0(set.ident, "_v2", set.y.max = 25))
paired.list = lapply(order.celltypes, function(i) Paired_t.test(celltype.x = i, prop.df = prop.df_idents))

# Create dataframe of all celltypes and their p-values after Mann-Whitney U-test
wilcox.celltypes.df = bind_rows(wilcox.list)
wilcox.celltypes.df$p.value_significant = wilcox.celltypes.df$array.wilcox.x. <= 0.05
wilcox.celltypes.df$p.adj_significant = wilcox.celltypes.df$array.wilcox.x.padj <= 0.05

paired.celltypes.df = bind_rows(paired.list)
paired.celltypes.df$p.value_significant = paired.celltypes.df$array.wilcox.x. <= 0.05
paired.celltypes.df$p.adj_significant = wilcox.celltypes.df$array.wilcox.x.padj <= 0.05


