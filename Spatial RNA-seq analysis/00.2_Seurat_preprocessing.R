#!/usr/bin/env Rscript

# Loading required packages
library(Seurat)
library(SeuratDisk)
library(anndata)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# QC cutoffs:
set_mt.cutoff = 20 
set_hb.cutoff = 1
set_nFeature.max = 10000 # Leave this higher
set_nFeature.min = 250
set_nCell_features = 3
set_ribo.cutoff = 5

# Flag to only run conversion
only.conversion = FALSE
do.plotting = TRUE

# Setting the input directory
Input.dir = "Data/0_Preprocessing_unfiltered_bin30_seurat/"

# Setting the output directory
if (dir.exists(path = "Output/0_Preprocessing") == FALSE) {
  print("Output/0_Preprocessing")
  dir.create(path = "Output/0_Preprocessing", recursive = TRUE)
  Output.dir = "Output/0_Preprocessing/"
} else if (dir.exists(path = "Output/0_Preprocessing") == TRUE) {
  print("Directory exists")
  Output.dir = "Output/0_Preprocessing/"
} else {
  print("Error with output directory")
}

# Generate a list of Stereoseq input directories
Input.list = list.files(Input.dir, recursive = FALSE, full.names = TRUE)

print("Running preprocessing on:")
for (name in Input.list) {
  print(name)
}

# Function to perform preprocessing and quality control
Preprocessing_QC <- function(seurat.x, mt.cutoff = set_mt.cutoff, ribo.cutoff = set_ribo.cutoff, hb.cutoff = set_hb.cutoff, 
                             nFeature.min = set_nFeature.min, nFeature.max = set_nFeature.max, nCell_feature_cutoff = set_nCell_features,
                             conversion.flag = only.conversion) {
  
  # Extract the samples name
  seurat.sample = sub(".*/([^/]+)$", "\\1", seurat.x)
  seurat.sample = sub("^(.*?)_", "", seurat.sample)
  seurat.sample = sub("_unfiltered_Seurat.h5ad", "", seurat.sample)
  print(paste0("Running preprocessing on ", seurat.sample))
  
  # Set the output directory for the sample
  if (dir.exists(path = paste0(Output.dir, seurat.sample)) == FALSE) {
    print(paste0(Output.dir, seurat.sample))
    dir.create(path = paste0(Output.dir, seurat.sample), recursive = TRUE)
    Output.dir.sample = paste0(Output.dir, seurat.sample, "/")
  } else if (dir.exists(path = paste0(Output.dir, seurat.sample)) == TRUE) {
    print("Directory exists")
    Output.dir.sample = paste0(Output.dir, seurat.sample, "/")
  } else {
    print("Error with output directory")
  }
  
  # Convert the h5ad to seurat
  # If only conversion is done, the script will exit the function
  check.Seurat = list.files(Output.dir.sample, pattern = "*.h5Seurat", full.names = TRUE)

  if (length(check.Seurat) == 1) {
    print("Seurat file exist, skipping conversion")
    seurat.x = LoadH5Seurat(check.Seurat, meta.data = FALSE)
  } else if (length(check.Seurat) == 0) {
    print("No Seurat found, run conversion")
    if (conversion.flag == TRUE) {
      print(paste("ONLY CONVERSION ON", seurat.sample))
      convert = h5ad_seurat_conversion(h5ad.x = seurat.x, seurat.output = Output.dir.sample, return.seurat = FALSE)
      print(paste("CONVERSION FINISHED OF", seurat.sample))
      return()
    } else if (conversion.flag == FALSE) {
      seurat.x = h5ad_seurat_conversion(h5ad.x = seurat.x, seurat.output = Output.dir.sample)
    } 
  } else if (length(check.Seurat) > 1) {
    print("Multiple Seurat files found in output. Resolve before continuing. Exiting preprocessing")
    return()
  }
  
  # Normalise and scale the data 
  DefaultAssay(seurat.x) <- "RNA"
  seurat.x = NormalizeData(seurat.x, normalization.method = "LogNormalize")
  seurat.x = FindVariableFeatures(seurat.x, selection.method = "vst", verbose = T)
  top10_var <- head(VariableFeatures(seurat.x), 10)
  all.genes = rownames(seurat.x)
  
  # Add the sample ident to the seurat object
  seurat.x$orig.ident = seurat.sample
  Idents(seurat.x) <- "orig.ident"
  DefaultAssay(seurat.x) <- "RNA"
  
  # Calculate mt% and hb%
  seurat.x[["percent.mt"]] <- PercentageFeatureSet(seurat.x, pattern = "^MT-", assay = "RNA")
  mt.genes <- rownames(seurat.x)[grep("^MT-",rownames(seurat.x))]
  
  seurat.x[["percent.hb"]] <- PercentageFeatureSet(seurat.x, pattern = "^HB[^(P)]", assay = "RNA")
  hb.genes <- rownames(seurat.x)[grep("^HB[^(P)]",rownames(seurat.x))]
  
  # Plotting prefiltering
  Idents(seurat.x) <- "orig.ident"
  print("Plot QC prefiltering")
  Plotting_QC = QC_plotting(x.QC = seurat.x, Output.dir = Output.dir.sample, stage = "prefiltering")
  
  # Perform filtering with set thresholds
  DefaultAssay(seurat.x) <- "RNA"
  
  # Select cells with the threshold
  selected_cells = WhichCells(seurat.x, expression = nFeature_RNA > nFeature.min &
                                nFeature_RNA < nFeature.max &
                                percent.mt < mt.cutoff &
                                percent.hb < hb.cutoff)
  selected_features = rownames(seurat.x) [Matrix::rowSums(seurat.x) > nCell_feature_cutoff]
  
  # Subsetting the cells based on selected cells
  print("Dimension of objects before filtering")
  print(dim(seurat.x))
  seurat.x = subset(seurat.x, features = selected_features, cells = selected_cells)
  print("Dimension of objects after filtering")
  print(dim(seurat.x))
  
  # Re-normalise and scale the data after filtering
  DefaultAssay(seurat.x) <- "RNA"
  seurat.x = NormalizeData(seurat.x, normalization.method = "LogNormalize")
  seurat.x = FindVariableFeatures(seurat.x, selection.method = "vst", verbose = T)
  top10_var <- head(VariableFeatures(seurat.x), 10)
  all.genes = rownames(seurat.x)
  seurat.x = ScaleData(seurat.x, verbose = T, features = all.genes)
  
  # QC plotting after filtering
  Plotting_QC = QC_plotting(x.QC = seurat.x, Output.dir = Output.dir.sample, stage = "filtered")
  
  # Perform scTransformation normalisation
  DefaultAssay(seurat.x) <- "RNA"
  seurat.x = SCTransform(seurat.x, verbose = T)
  
  # Do PCA on scTranformed object
  seurat.x = RunPCA(seurat.x)
  
  # Plotting PCA results
  PCA_plot = DimPlot(seurat.x, reduction = "pca")
  ggsave2(plot = PCA_plot, filename = paste0(Output.dir.sample, "scTransform_PCA_plot.pdf"))
  
  PCA_heatmap = DimHeatmap(seurat.x, dims = 1:15, cells = 500, balanced = TRUE)
  ggsave2(plot = PCA_heatmap, filename = paste0(Output.dir.sample, "scTransform_PCA_heatmap.pdf"))
  
  # Determine he dimensionality of the data
  elbow_plot = ElbowPlot(seurat.x, reduction = "pca", ndims = 25)
  ggsave2(plot = elbow_plot, filename = paste0(Output.dir.sample, "scTransform_elbow_plot.pdf"))
  
  # Do UMAP on seurat object
  seurat.x = RunUMAP(seurat.x, dims = 1:30)
  
  # Run clustering on the object
  seurat.x = FindNeighbors(object = seurat.x, reduction = "pca", dims = 1:30, verbose = FALSE)
  seurat.x = FindClusters(object = seurat.x, resolution = 0.7, verbose = FALSE)
  
  # Plotting UMAP
  UMAP_plot = DimPlot(seurat.x, label = TRUE, repel = TRUE)
  ggsave2(plot = UMAP_plot, filename = paste0(Output.dir.sample, "scTransform_UMAP.pdf"))
  
  # Plotting UMAP on spatial
  UMAP_plot = DimPlot(seurat.x, label = FALSE, repel = FALSE, reduction = "spatial")
  ggsave2(plot = UMAP_plot, filename = paste0(Output.dir.sample, "scTransform_UMAP_spatial.pdf"))
  
  # Reset deafault to RNA
  DefaultAssay(seurat.x) <- "RNA"
  
  # Compute cell cycle scoring
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  seurat.x = CellCycleScoring(seurat.x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  # Plot the cell cycle score
  # Visualize the distribution of cell cycle markers across
  VlnPlot(seurat.x, features = c("S.Score","G2M.Score"), group.by = "seurat_clusters", pt.size = 0)
  ggsave2(paste0(Output.dir.sample, "Cell.Cycle_Phase_VlnPlot_filtered.pdf"))
  
  # Plot cell cycle phases
  Phase_UMAP = DimPlot(seurat.x, group.by = "Phase")
  ggsave2(plot = Phase_UMAP, filename = paste0(Output.dir.sample, "Cell-Cycle_Phase_UMAP_filtered.pdf"))
  
  # Plotting additional UMAPs before saving
  FeaturePlot(seurat.x, features = "percent.mt")
  ggsave2(paste0(Output.dir.sample, "UMAP_percent_mt_filtered.pdf"))
  
  FeaturePlot(seurat.x, features = "percent.hb")
  ggsave2(paste0(Output.dir.sample, "UMAP_percent_hb_filtered.pdf"))
  
  FeaturePlot(seurat.x, features = "nFeature_SCT")
  ggsave2(paste0(Output.dir.sample, "UMAP_nFeature_SCT_filtered.pdf"))
  
  FeaturePlot(seurat.x, features = "nFeature_RNA")
  ggsave2(paste0(Output.dir.sample, "UMAP_nFeature_RNA_filtered.pdf"))
  
  FeaturePlot(seurat.x, features = "nCount_SCT")
  ggsave2(paste0(Output.dir.sample, "UMAP_nCount_SCT_filtered.pdf"))
  
  FeaturePlot(seurat.x, features = "nCount_RNA")
  ggsave2(paste0(Output.dir.sample, "UMAP_nCount_RNA_filtered.pdf"))
  
  FeaturePlot(seurat.x, features = "S.Score")
  ggsave2(paste0(Output.dir.sample, "UMAP_S_score_filtered.pdf"))
  
  FeaturePlot(seurat.x, features = "G2M.Score")
  ggsave2(paste0(Output.dir.sample, "UMAP_G2M_Score_filtered.pdf"))

  # Save the filtered Seurat object
  SaveH5Seurat(object = seurat.x, filename = paste0(Output.dir.sample, seurat.sample, ".h5Seurat"), overwrite = TRUE)
  
  # Remove to clear memory
  seurat.x = NULL
  #return(x.seurat)
  
}

#### Function block ####
h5ad_seurat_conversion = function(h5ad.x = seurat.x, seurat.output = Output.dir.sample, return.seurat = TRUE) {
  
  # Extract the samples name
  seurat.sample = sub(".*/([^/]+)$", "\\1", h5ad.x)
  seurat.sample = sub("^(.*?)_", "", seurat.sample)
  seurat.sample = sub("_unfiltered_Seurat.h5ad", "", seurat.sample)
  print(paste0("Running preprocessing on ", seurat.sample))
  
  # Set the output directory for the sample
  if (dir.exists(path = paste0(Output.dir, seurat.sample)) == FALSE) {
    print(paste0(Output.dir, seurat.sample))
    dir.create(path = paste0(Output.dir, seurat.sample), recursive = TRUE)
    Output.dir.sample = paste0(Output.dir, seurat.sample, "/")
  } else if (dir.exists(path = paste0(Output.dir, seurat.sample)) == TRUE) {
    print("Directory exists")
    Output.dir.sample = paste0(Output.dir, seurat.sample, "/")
  } else {
    print("Error with output directory")
  }
  
  # Convert the sample from h5ad to h5Seurat
  seurat.h5Seurat = paste0(Output.dir.sample, seurat.sample, ".h5Seurat")
  Convert(source = h5ad.x, dest = seurat.h5Seurat, overwrite = TRUE)
  h5ad.x = NULL
  
  # Load the newly generated h5Seurat object with meta.data
  seurat.x = list.files(Output.dir.sample, pattern = "*.h5Seurat", full.names = TRUE)
  print(paste("Loading converted h5Seurat", seurat.x))
  seurat.x = LoadH5Seurat(seurat.x, meta.data=FALSE) # meta.data is false to be able to be loaded
  print("Conversion has finished")
  
  if (return.seurat == TRUE) {
    return(seurat.x)
  } 
  
}



QC_plotting <- function(x.QC = seurat.x, stage, Output.dir = Output.dir.sample, FeaturePlot_col = c("#FFFFB2", "#B10026"), 
                        mt.cutoff = set_mt.cutoff, 
                        hb.cutoff = set_hb.cutoff, 
                        nFeature.min = set_nFeature.min, 
                        nFeature.max = set_nFeature.max) {
  
  print("Plotting QC")
  
  # Visualize QC metrics as a violin plot
  VlnPlot(x.QC, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0) +
    theme_cowplot()
  ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,".pdf"))
  
  if (nFeature.min == 0 && mt.cutoff == 0) {
    
    # Violinplots
    VlnPlot(x.QC, features = "nFeature_RNA", pt.size = 0) + 
      theme_cowplot()
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_nFeatures.pdf"))
    
    VlnPlot(x.QC, features = "percent.mt", pt.size = 0) + 
      theme_cowplot()
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_percent-mt.pdf"))
    
    VlnPlot(x.QC, features = "percent.hb", pt.size = 0) + 
      theme_cowplot()
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_percent-hb.pdf"))
    
    # FeatureScatter plot nCount_RNA vs. percent.mt & nFeature_RNA
    plot_nFeature = FeatureScatter(x.QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      theme_cowplot()
    ggsave2(plot = plot_nFeature, paste0(Output.dir, "Scatterplot_QC_nFeature_", stage,".pdf"))
    plot_mt = FeatureScatter(x.QC, feature1 = "nCount_RNA", feature2 = "percent.mt")  +
      theme_cowplot()
    ggsave2(plot = plot_mt, paste0(Output.dir, "Scatterplot_QC_mtDNA_", stage,".pdf"))
  } else if (nFeature.min > 0 && mt.cutoff > 0) {
    
    # Violinplots
    VlnPlot(x.QC, features = "nFeature_RNA", pt.size = 0) + 
      geom_hline(aes(yintercept = nFeature.min, linetype = "dashed")) +
      theme_cowplot()
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_nFeatures.pdf"))
    
    VlnPlot(x.QC, features = "percent.mt", pt.size = 0) + 
      geom_hline(aes(yintercept = mt.cutoff), linetype = "dashed") +
      theme_cowplot()
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_percent-mt.pdf"))
    
    VlnPlot(x.QC, features = "percent.hb", pt.size = 0) + 
      geom_hline(aes(yintercept = hb.cutoff), linetype = "dashed") +
      theme_cowplot()
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_percent-hb.pdf"))
    
    # FeatureScatter plot nCount_RNA vs. percent.mt & nFeature_RNA
    plot_nFeature = FeatureScatter(x.QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      geom_hline(aes(yintercept = nFeature.min, linetype = "dashed")) +
      theme_cowplot()
    ggsave2(plot = plot_nFeature, paste0(Output.dir, "Scatterplot_QC_nFeature_", stage,".pdf"))
    plot_mt = FeatureScatter(x.QC, feature1 = "nCount_RNA", feature2 = "percent.mt")  +
      geom_hline(aes(yintercept = mt.cutoff, linetype = "dashed")) +
      theme_cowplot()
    ggsave2(plot = plot_mt, paste0(Output.dir, "Scatterplot_QC_mtDNA_", stage,".pdf"))
  }
  
  # Generate featureplots
  spatial_nFeature = FeaturePlot(x.QC, features = "nFeature_RNA", reduction = "spatial", cols = FeaturePlot_col)
  ggsave2(plot = spatial_nFeature, paste0(Output.dir, "FeaturePlot_nFeature_", stage,".pdf"))
  spatial_nCount = FeaturePlot(x.QC, features = "nCount_RNA", reduction = "spatial", cols = FeaturePlot_col)
  ggsave2(plot = spatial_nFeature, paste0(Output.dir, "FeaturePlot_nCount_", stage,".pdf"))
  spatial_mtDNA = FeaturePlot(x.QC, features = "percent.mt", reduction = "spatial", cols = FeaturePlot_col)
  ggsave2(plot = spatial_mtDNA, paste0(Output.dir, "FeaturePlot_mtDNA_", stage,".pdf"))
  
  x.QC = NULL
  
}

# Run pre-processing on all samples. 
Do_preprocessing = lapply(Input.list, Preprocessing_QC)

