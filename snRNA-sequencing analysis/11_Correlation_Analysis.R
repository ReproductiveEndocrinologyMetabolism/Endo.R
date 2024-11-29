#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(reshape2)
library(dplyr)
library(tidyr)

# Set variables
corr.cutoff = 0.5
p.cutoff = 0.05
Group.subset = c("Control", "PCOS_W0", "PCOS") # In Clin it is PCOS, in Seurat it is PCOS_W0
Clin.cols = c("FerrimanGallwey", 
              "DHEA_pgml", 
              "Adion_pgml",
              "T_pgml",
              "E2_pgml",
              "Prog_pgml",
              "SHBG_nmol_L",
              "FAI",
              "LH", 
              "FSH",
              "HomaIR")
Clin.cols = c("FerrimanGallwey", 
              "DHEA_pgml", 
              "Adion_pgml",
              "LH", 
              "HomaIR")

# Set to load averaged seruat or not
load.average = FALSE
load.seurat = TRUE
save.average = FALSE
#Seurat.path = "Data/CPM_Averaged_Seurat/" # Averaged object
Normalization.method = "Average" # CPM or LogNormalize or Average
Seurat.path = "Data/Normalized_Seurat_Subset/" # Seurat object

# Set the path to clinical variables table and common DEG table
Clin.df = read.xlsx(xlsxFile = "Data/Clinical_data/Clinical_data.xlsx") # Clinical data
DEG.df = read.xlsx(xlsxFile = "Data/DEG_tables/DEG_table_all_celltypes_Control_PCOS.xlsx") # DEG list combined cell types

# List the subset cell types
Celltypes.list = list.files(path = Seurat.path)
#Celltypes.list = Celltypes.list[-28] # Remove undefined celltype cluster

# Generate or set the output directory
if (dir.exists(path = paste0("Output/12_Correlation_Analysis_", Normalization.method)) == FALSE) {
  print(paste0("Generating output directory Output/12_Correlation_Analysis"))
  dir.create(path = paste0("Output/12_Correlation_Analysis_", Normalization.method), recursive = TRUE)
  Output.dir = paste0("Output/12_Correlation_Analysis_", Normalization.method, "/")
} else if (dir.exists(path = paste0("Output/12_Correlation_Analysis_", Normalization.method)) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/12_Correlation_Analysis_", Normalization.method, "/")
} else {
  print("Error with output directory")
}

# Extract target groups from clinical data
Clin.df = Clin.df[Clin.df$Group %in% Group.subset,]

# Check if the clinical variables are normally distributed with a shapiro-wilks test.
# Do the test group wise and together
Clin_cols.shapiro = which(colnames(Clin.df) %in% Clin.cols)

# Initialize result dataframe
Shapiro.df = data.frame(
  Variable = colnames(Clin.df)[Clin_cols.shapiro],
  Control = numeric(length(Clin_cols.shapiro)),
  PCOS = numeric(length(Clin_cols.shapiro)),
  Control_PCOS = numeric(length(Clin_cols.shapiro))
)

# Run Shapiro-Wilk test and collect p-values
for (i in seq_along(Clin_cols.shapiro)) {
  col <- Clin_cols.shapiro[i]
  variable <- colnames(Clin.df)[col]
  
  control_values <- Clin.df %>% filter(Group == "Control") %>% pull(col)
  pcos_values <- Clin.df %>% filter(Group == "PCOS") %>% pull(col)
  control_pcos_values <- Clin.df %>% pull(col)
  
  Shapiro.df$Control[i] <- shapiro.test(control_values)$p.value
  Shapiro.df$PCOS[i] <- shapiro.test(pcos_values)$p.value
  Shapiro.df$Control_PCOS[i] <- shapiro.test(control_pcos_values)$p.value
}

write.xlsx(x = Shapiro.df, file = paste0(Output.dir, "Shapiro_Wilks_table.xlsx"))

# Generate output list containing all correlation for all, Control and PCOS
output.list_All_Corr_Control_PCOS = list()
output.list_All_Corr_Control = list()
output.list_All_Corr_PCOS = list()

# Generate output list containing correlations passing the cutoff
output.list_Cutoff_Corr_Control_PCOS = list()
output.list_Cutoff_Corr_Control = list()
output.list_Cutoff_Corr_PCOS = list()

# Generate output list for correlations passing the cutoff found in both groups
output.list_Cutoff_Corr_Overlap = list()

# Loop over every cell type
for (celltype.x in Celltypes.list) {
  
  if (load.seurat == TRUE && load.average == FALSE) {
    
    # Extract the cell type
    CT.name = gsub("Endo_All_(.*)_subset_seurat.h5seurat", "\\1", celltype.x)
    print(CT.name)
    
    # Subset the main DEG table for current cell type
    DEG.ct = DEG.df[DEG.df$Subtype == CT.name, ]
    
    # If there are less than 3 DEG, skip to next celltype
    if (nrow(DEG.ct) < 3) {
      print("Too few DEGs")
      next
    }
    
    # Load the subset cell type seurat object
    print("Loading the Seurat object")
    seurat.ct = LoadH5Seurat(file = paste0(Seurat.path, celltype.x), assays = "RNA")
    
    # Subsetting for selected group
    seurat.ct = subset(x = seurat.ct, subset = Group_Stage %in% Group.subset)
    
    # Set ident to orig.ident before running average expression.
    # Gives you average expression per sample
    Idents(seurat.ct) = "orig.ident"
    
    # Test aggregating data to pseudobulk with CPM
    if (Normalization.method == "CPM") {
      seurat.ct = AggregateExpression(object = seurat.ct, assays = "RNA", normalization.method = "RC", 
                                      scale.factor = 1e6, return.seurat = TRUE)
    } else if (Normalization.method == "LogNormalize") {
      seurat.ct = AggregateExpression(object = seurat.ct, assays = "RNA", normalization.method = "LogNormalize", 
                                      scale.factor = 10000, return.seurat = TRUE)
    } else if (Normalization.method == "Average") {
      seurat.ct = AverageExpression(object = seurat.ct, assays = "RNA", layer = "data", return.seurat = TRUE)
    }
    
    # Save the averaged object
    if (save.average == TRUE) {
      SaveH5Seurat(object = seurat.ct, filename = paste0(Output.dir, "Endo_Baseline_", CT.name, "_", Normalization.method, "_average.h5seurat"), overwrite = TRUE)
    }
    
    # Convert to a dataframe
    seurat.ct = data.frame(seurat.ct@assays$RNA$data)

  } else if (load.seurat == FALSE && load.average == TRUE) { # EDIT THIS LATER
    
    # Extract the cell type
    CT.name = gsub(paste0("Endo_All_(.*)_", Normalization.method, "_average.h5seurat"), "\\1", celltype.x)
    print(CT.name)
    
    # Subset the main DEG table for current cell type
    DEG.ct = DEG.df[DEG.df$Subtype == CT.name, ]
    
    # If there are less than 3 DEG, skip to next celltype
    if (nrow(DEG.ct) < 3) {
      print("Too few DEGs")
      next
    }
    
    # Load the averaged seurat object
    print("Loading the averaged object")
    seurat.ct = LoadH5Seurat(file = paste0(Seurat.path, celltype.x))
    
    # Subsetting for selected group
    seurat.ct = subset(x = seurat.ct, subset = Group_Stage %in% Group.subset)
    
    # Extract assay data
    df.count = data.frame(seurat.ct@assays$RNA$counts)
    df.data = data.frame(seurat.ct@assays$RNA$data)
    df.scale= data.frame(seurat.ct@assays$RNA$scale.data)
    sum(colSums(df.count))
    colSums(df.data)
    colSums(df.scale)
    
    df.alt.count = data.frame(seurat.ct@assays$RNA@counts)
    colSums(df.alt.count)
    
    if (Normalization.method == "CPM") {
      seurat.ct = data.frame(seurat.ct@assays$RNA$data)
    } else if (Normalization.method == "Log1p") {
      seurat.ct = data.frame(seurat.ct@assays$RNA$counts)
    }
    
    
  }
  
  # Subset the average dataframe for only DEG's
  avg.df = seurat.ct[rownames(seurat.ct) %in% DEG.ct$gene,]
  
  # Rename colnames of average df to fit clinical dataframe
  colnames(avg.df) <- sub("Ctrl_KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- sub("LS_KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- sub("Met_KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- gsub("\\_W[0-9]+", "", colnames(avg.df))
  
  colnames(avg.df) <- sub("Ctrl.KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- sub("LS.KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- sub("Met.KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- gsub("\\.W[0-9]+", "", colnames(avg.df))
  
  # For-loop to run linear regression
  DEG.vec = rownames(avg.df)
  
  for (gene.x in DEG.vec) {
    
    # Extract the selected DEG
    gene.DEG = avg.df[rownames(avg.df) == gene.x,]
    
    # Convert to long format
    gene.DEG = melt(gene.DEG, variable.name = "ID", value.name = Normalization.method)
    
    # Merge the DEG and clinical var data frames based on id
    Clin_DEG.df <- merge(Clin.df, gene.DEG, by = "ID", all.x = TRUE)
    
    # Subset columns to correlate
    Corr_var.df_all = Clin_DEG.df[,c(Clin.cols, Normalization.method)] # Edit this depending on your Clinical table
    Corr_var.df_Control = Clin_DEG.df[Clin_DEG.df$Group == "Control", c(Clin.cols, Normalization.method)]
    Corr_var.df_PCOS = Clin_DEG.df[Clin_DEG.df$Group == "PCOS", c(Clin.cols, Normalization.method)]
    
    # Calculate the correlations
    Corr_var.df_all <- cor(Corr_var.df_all, use = "complete.obs", method = "spearman")
    Corr_var.df_Control <- cor(Corr_var.df_Control, use = "complete.obs", method = "spearman")
    Corr_var.df_PCOS <- cor(Corr_var.df_PCOS, use = "complete.obs", method = "spearman")
    
    # Extract correlations with gene expression
    Corr.expression_all <-  Corr_var.df_all[Normalization.method,]
    Corr.expression_Control <-  Corr_var.df_Control[Normalization.method,]
    Corr.expression_PCOS <-  Corr_var.df_PCOS[Normalization.method,]
    
    # Generate a dataframe of the output, all correlations
    output_all.df = data.frame(Celltype = CT.name, "Gene" = gene.x, t(as.data.frame(Corr.expression_all)))
    output_Control.df = data.frame(Celltype = CT.name, "Gene" = gene.x, t(as.data.frame(Corr.expression_Control)))
    output_PCOS.df = data.frame(Celltype = CT.name, "Gene" = gene.x, t(as.data.frame(Corr.expression_PCOS)))
    
    # Append the dataframes to output lists, all correlations
    output.list_All_Corr_Control_PCOS[[paste0(CT.name, "_", gene.x)]] = output_all.df
    output.list_All_Corr_Control[[paste0(CT.name, "_", gene.x)]] = output_Control.df
    output.list_All_Corr_PCOS[[paste0(CT.name, "_", gene.x)]] = output_PCOS.df
    
    # Extract correlation over cutoff
    Corr.expression_all = Corr.expression_all[Corr.expression_all >= corr.cutoff | Corr.expression_all <=  -corr.cutoff]
    Corr.expression_all = Corr.expression_all[names(Corr.expression_all) != Normalization.method]
    Corr.expression_all = Corr.expression_all[!is.na(Corr.expression_all)]
    
    Corr.expression_Control = Corr.expression_Control[Corr.expression_Control >= corr.cutoff | Corr.expression_Control <=  -corr.cutoff]
    Corr.expression_Control = Corr.expression_Control[names(Corr.expression_Control) != Normalization.method]
    Corr.expression_Control = Corr.expression_Control[!is.na(Corr.expression_Control)]
    
    Corr.expression_PCOS = Corr.expression_PCOS[Corr.expression_PCOS >= corr.cutoff | Corr.expression_PCOS <=  -corr.cutoff]
    Corr.expression_PCOS = Corr.expression_PCOS[names(Corr.expression_PCOS) != Normalization.method]
    Corr.expression_PCOS = Corr.expression_PCOS[!is.na(Corr.expression_PCOS)]
    
    # Check if the there are variables passing the cutoff. Then generate output
    if (length(Corr.expression_all) > 0) {
      
      # Generate the dataframe and append to output list
      output.df = data.frame(Celltype = CT.name, "Gene" = gene.x, t(as.data.frame(Corr.expression_all)))
      output.list_Cutoff_Corr_Control_PCOS[[paste0(CT.name, "_", gene.x)]] = output.df
    }
    if (length(Corr.expression_Control) > 0) {
      
      # Generate the dataframe and append to output list
      output.df = data.frame(Celltype = CT.name, "Gene" = gene.x, t(as.data.frame(Corr.expression_Control)))
      output.list_Cutoff_Corr_Control[[paste0(CT.name, "_", gene.x)]] = output.df
    }
    if (length(Corr.expression_PCOS) > 0) {
      
      # Generate the dataframe and append to output list
      output.df = data.frame(Celltype = CT.name, "Gene" = gene.x, t(as.data.frame(Corr.expression_PCOS)))
      output.list_Cutoff_Corr_PCOS[[paste0(CT.name, "_", gene.x)]] = output.df
    }
    
    # Check if and which variables correlations overlap between Control and PCOS
    if (length(Corr.expression_Control) > 0 & length(Corr.expression_PCOS) > 0) {
      
      # Check for overlapping variables
      Corr.expression_Overlap = intersect(names(Corr.expression_Control), names(Corr.expression_PCOS))
      
      if (length(Corr.expression_Overlap) > 0) {
        
        Corr.expression_all_overlap <-  Corr_var.df_all[Normalization.method,]
        Corr.expression_all_overlap = Corr.expression_all_overlap[Corr.expression_all_overlap >= corr.cutoff | Corr.expression_all_overlap <=  -corr.cutoff]
        
        # Extract the correlation for all correlations
        Corr.expression_all_overlap = Corr.expression_all_overlap[names(Corr.expression_all_overlap) %in% Corr.expression_Overlap]
        
        # Generate a dataframe and append to output
        output.df = data.frame(Celltype = CT.name, "Gene" = gene.x, t(as.data.frame(Corr.expression_all_overlap)))
        output.list_Cutoff_Corr_Overlap[[paste0(CT.name, "_", gene.x)]] = output.df
        
      }
      
      
    }
    
    
  }
  
}

# Clean the memory
seurat.ct = NULL

# Combined the lists to dataframes, all correlations
output.df_All_Corr_Control_PCOS = bind_rows(output.list_All_Corr_Control_PCOS)
output.df_All_Corr_Control = bind_rows(output.list_All_Corr_Control)
output.df_All_Corr_PCOS = bind_rows(output.list_All_Corr_PCOS)

# Combined the lists to dataframes, correlations passing the cutoff
output.df_Cutoff_Corr_Control_PCOS = bind_rows(output.list_Cutoff_Corr_Control_PCOS)
output.df_Cutoff_Corr_Control = bind_rows(output.list_Cutoff_Corr_Control)
output.df_Cutoff_Corr_PCOS = bind_rows(output.list_Cutoff_Corr_PCOS)

# Combined the lists to dataframes, correlations passing the cutoff found in both groups
output.df_Cutoff_Corr_Overlap = bind_rows(output.list_Cutoff_Corr_Overlap)

# Save the output dataframes
write.xlsx(x = output.df_All_Corr_Control_PCOS, file = paste0(Output.dir, "Spearman_Corr_All_Control_PCOS.xlsx"))
write.xlsx(x = output.df_All_Corr_Control, file = paste0(Output.dir, "Spearman_Corr_All_Control.xlsx"))
write.xlsx(x = output.df_All_Corr_PCOS, file = paste0(Output.dir, "Spearman_Corr_All_PCOS.xlsx"))

write.xlsx(x = output.df_Cutoff_Corr_Control_PCOS, file = paste0(Output.dir, "Spearman_Corr_Cutoff_", corr.cutoff, "_Control_PCOS.xlsx"))
write.xlsx(x = output.df_Cutoff_Corr_Control, file = paste0(Output.dir, "Spearman_Corr_Cutoff_", corr.cutoff, "_Control.xlsx"))
write.xlsx(x = output.df_Cutoff_Corr_PCOS, file = paste0(Output.dir, "Spearman_Corr_Cutoff_", corr.cutoff, "_PCOS.xlsx"))

write.xlsx(x = output.df_Cutoff_Corr_Overlap, file = paste0(Output.dir, "Spearman_Corr_Cutoff_", corr.cutoff, "_Overlap.xlsx"))

# Pivot the overlap table for downstream correlation test. Remove NA
Overlapped.Corr <- output.df_Cutoff_Corr_Overlap %>%
  pivot_longer(cols = -c(Celltype, Gene), names_to = "Measurement", values_to = "Correlation")
Overlapped.Corr = na.omit(Overlapped.Corr)

# Filter for Control-PCOS correlation over cutoff
Overlapped.Corr = Overlapped.Corr[Overlapped.Corr$Correlation >= corr.cutoff | Overlapped.Corr$Correlation <=  -corr.cutoff,]

# Pivot the Control_PCOS table for downstream correlation test. Remove NA
Control_PCOS.Corr <- output.df_Cutoff_Corr_Control_PCOS %>%
  pivot_longer(cols = -c(Celltype, Gene), names_to = "Measurement", values_to = "Correlation")
Control_PCOS.Corr = na.omit(Control_PCOS.Corr)

# Filter for Control-PCOS correlation over cutoff
Control_PCOS.Corr = Control_PCOS.Corr[Control_PCOS.Corr$Correlation >= corr.cutoff | Control_PCOS.Corr$Correlation <=  -corr.cutoff,]

# Check which celltypes are in the downstream tables
Overlap.CT = unique(output.df_Cutoff_Corr_Overlap$Celltype)
Control_PCOS.CT = unique(output.df_Cutoff_Corr_Control_PCOS$Celltype)

if (load.seurat == TRUE) {
  CT.names <- sapply(Celltypes.list, function(x) gsub("Endo_All_(.*)_subset_seurat\\.h5seurat", "\\1", x))
  names(Celltypes.list) = CT.names
} else if (load.average == TRUE) {
  CT.names <- sapply(Celltypes.list, function(x) gsub(paste0("Endo_All_(.*)_", Normalization.method, "_average.h5seurat"), "\\1", x))
  names(Celltypes.list) = CT.names
}

##### DECIDE TO RUN DOWNSTREAM ON EITHER OVERLAPPED OR MERGED TABLES I.E
##### OVERLAPPED OR CONTROL-PCOS 
Select.CT = Control_PCOS.CT # Overlap.CT OR Control_PCOS.CT
#Overlapped.Corr.df = as.data.frame(Overlapped.Corr)
Corr.df = as.data.frame(Control_PCOS.Corr) # Overlapped.Corr OR Control_PCOS.Corr


# Filter for celltypes present in overlap dataframe
Celltypes.list = Celltypes.list[names(Celltypes.list) %in% Select.CT]

# Second loop to run correlation tests on overlapped variables to detect significant correlation
corr.output.list = list()
plot.output.list = list()

for (celltype.x in Celltypes.list) {
  
  if (load.seurat == TRUE && load.average == FALSE) {
    
    # Extract the cell type
    CT.name = gsub("Endo_All_(.*)_subset_seurat.h5seurat", "\\1", celltype.x)
    print(CT.name)
    
    # Load the subset cell type seurat object
    print("Loading the Seurat object")
    seurat.ct = LoadH5Seurat(file = paste0(Seurat.path, celltype.x), assays = "RNA")
    
    # Subsetting for selected group
    seurat.ct = subset(x = seurat.ct, subset = Group_Stage %in% Group.subset)
    
    # Set ident to orig.ident before running average expression.
    # Gives you average expression per sample
    Idents(seurat.ct) = "orig.ident"
    
    # Do average expression
    if (Normalization.method == "CPM") {
      seurat.ct = AggregateExpression(object = seurat.ct, assays = "RNA", normalization.method = "RC", 
                                      scale.factor = 1e6, return.seurat = TRUE)
    } else if (Normalization.method == "LogNormalize") {
      seurat.ct = AggregateExpression(object = seurat.ct, assays = "RNA", normalization.method = "LogNormalize", 
                                      scale.factor = 10000, return.seurat = TRUE)
    } else if (Normalization.method == "Average") {
      seurat.ct = AverageExpression(object = seurat.ct, assays = "RNA", layer = "data", return.seurat = TRUE)
    }
    avg.df = data.frame(seurat.ct@assays$RNA$data)
    
  } else if (load.seurat == FALSE && load.average == TRUE) { # EDIT THIS LATER
    
    # Extract the cell type
    CT.name = gsub(paste0("Endo_All_(.*)_", Normalization.method, "_average.h5seurat"), "\\1", celltype.x)
    print(CT.name)
    
    # Load the averaged seurat object
    print("Loading the averaged object")
    avg.df = LoadH5Seurat(file = paste0(Seurat.path, celltype.x))
    avg.df = data.frame(avg.df@assays$RNA$data)
    
  }
  
  # Subset the main DEG table for current cell type
  DEG.ct = DEG.df[DEG.df$Subtype == CT.name, ]
  
  # Rename colnames of average df to fit clinical dataframe
  colnames(avg.df) <- sub("Ctrl_KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- sub("LS_KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- sub("Met_KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- gsub("\\_W[0-9]+", "", colnames(avg.df))
  
  colnames(avg.df) <- sub("Ctrl.KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- sub("LS.KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- sub("Met.KS", "ks",  colnames(avg.df))
  colnames(avg.df) <- gsub("\\.W[0-9]+", "", colnames(avg.df))
  
  # Loop over only DEG's found in overlapped dataframe
  Corr.CT = Corr.df[Corr.df$Celltype == CT.name,]
  DEG.vec = unique(Corr.CT$Gene)
  
  for (gene.x in DEG.vec) {
    
    print(paste("Correlation on", CT.name, "and", gene.x))
    
    # Extract the selected DEG
    gene.DEG = avg.df[rownames(avg.df) == gene.x,]
    
    # Convert to long format
    gene.DEG = melt(gene.DEG, variable.name = "ID", value.name = Normalization.method)
    
    # Merge the DEG and clinical var data frames based on id
    Clin_DEG.df <- merge(Clin.df, gene.DEG, by = "ID", all.x = TRUE)
    
    # Subset the celltype specific overlapped dataframe for the current gene
    Corr.gene = Corr.CT[Corr.CT$Gene == gene.x, ]
    
    # Loop over overlapped correlated values in selected celltype and gene
    Corr.vars = Corr.gene$Measurement
    
    #Corr.x = "HomaIR"
    for (Corr.x in Corr.vars) {
      
      # Perform the correlation test
      cor.out = cor.test(Clin_DEG.df[[Normalization.method]], Clin_DEG.df[[Corr.x]], method = "spearman")
      cor.rho = cor.out$estimate
      cor.p = cor.out$p.value
      
      output.df = data.frame(Celltype = CT.name, "Gene" = gene.x, Measurment = Corr.x, 
                             Correlation_Rho = cor.rho, P.value = cor.p)
      corr.output.list[[paste0(CT.name, "_", gene.x, "_", Corr.x)]] = output.df
      
      if (cor.p <= p.cutoff) {
        plot.x = ggscatter(Clin_DEG.df, x = Corr.x, y = Normalization.method, 
                           color = "Group", shape = "Group", size = 10, # Old size is 4
                           palette = c("#A0A0A0", "#D098B2"),
                           add = "reg.line",  # Add regressin line
                           conf.int = FALSE,
                           add.params = list(color = "red"), # Customize reg. line
                           cor.coef = TRUE, 
                           cor.method = "spearman",
                           xlab = Corr.x, ylab = paste(gene.x, Normalization.method))
        plot.x
        plot.output.list[[paste0(CT.name, "_", gene.x, "_", Corr.x)]] = plot.x
      }
      
    }
    
  }
  
}

# Combined the list to dataframe, correlation tested 
corr.output.df = bind_rows(corr.output.list)

# Save the output dataframes
write.xlsx(x = corr.output.df, file = paste0(Output.dir, "Spearman_Corr_test_Baseline.xlsx"))

# Save all the plots in the Plots folder
dir.create(path = paste0(Output.dir, "Plots"), recursive = TRUE)
Output_plot.dir = paste0(Output.dir, "Plots/")

# Function to save plots
save_plot <- function(plot_name, plot) {
  file_name <- paste0(Output_plot.dir, plot_name, ".pdf")
  ggsave(filename = file_name, plot = plot, device = "pdf", dpi = 300)
}

# Apply the function to each plot in the list
lapply(names(plot.output.list), function(plot_name) {
  save_plot(plot_name, plot.output.list[[plot_name]])
})
