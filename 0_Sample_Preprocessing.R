library(Seurat)
library(remotes)
library(cowplot)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(hdf5r)
library(SeuratDisk) 
library(scDblFinder)
library(future)

#Install latest scTransform
#devtools::install_github("satijalab/sctransform", ref = "develop")

 # Setting the input directory
Input_dir = "Data/CellRanger_Count/"
Input.list = list.files(Input_dir)

# QC cutoffs:
mt.cutoff = 10 
ribo.cutoff = 5
hb.cutoff = 1
nFeature.max = 5000 # Leave this higher
nFeature.min = 500
nCell_features = 3
log2_norm_flag = FALSE
scT_norm_flag = TRUE

pre_processing.fun<- function(x, mt.cutoff = 10, ribo.cutoff = 5, 
                              hb.cutoff = 2, nFeature.max = 15000, 
                              nFeature.min = 500, nCell_features = 3) {
  
  #Generate output directory
  if (dir.exists(path = paste0("Output/0_Preprocessed/", x)) == FALSE) {
    print(paste0("Generating output directory ", "Output/", x))
    dir.create(path = paste0("Output/0_Preprocessed/", x), recursive = TRUE)
    Output.dir = paste0("Output/0_Preprocessed/", x, "/")
  } else if (dir.exists(path = paste0("Output/0_Preprocessed/", x)) == TRUE) {
    print("Directory exists")
    Output.dir = paste0("Output/0_Preprocessed/", x, "/")
  } else {
    print("Error with output directory")
  }
  
  # Generating the Seurat object
  Input.x = paste0(Input_dir, x)
  
  if (length(list.files(paste0(Input_dir, x), pattern = "raw.h5seurat", full.names = TRUE)) == 1) {
    print(paste0("Loading ", list.files(paste0(Input_dir, x), pattern = "_raw.h5seurat", full.names = FALSE)))
    x.seurat = LoadH5Seurat(file = list.files(paste0(Input_dir, x), pattern = "_raw.h5seurat", full.names = TRUE))
  } else if (length(list.files(Input_dir, pattern = "filtered_feature_bc_matrix", full.names = TRUE)) == 1) {
    #Generating h5 seurat object
    print("Generating seurat object")
    x.data <- Read10X(data.dir = list.files(Input_dir, pattern = "filtered_feature_bc_matrix", full.names = TRUE))
    x.seurat <- CreateSeuratObject(counts = x.data, project = x, min.cells = 3, min.features = 200)
    # Save raw seurat object
    print(paste0("Saving as ", x, "_Endometrium_raw.h5seurat"))
    SaveH5Seurat(x.seurat, paste0(Output.dir, x, "_Endometrium_raw.h5seurat"), overwrite = TRUE)
  } else {
    print("Error in data directory")
  }
  
  ##### Quality control and selecting cells for further analysis ####
  x.seurat[["percent.mt"]] <- PercentageFeatureSet(x.seurat, pattern = "^MT-")
  mt.genes <- rownames(x.seurat)[grep("^MT-",rownames(x.seurat))]
  
  x.seurat[["percent.ribo"]] <- PercentageFeatureSet(x.seurat, pattern = "^RP[SL]")
  ribo.genes <- rownames(x.seurat)[grep("^RP[SL]",rownames(x.seurat))]
  
  x.seurat[["percent.hb"]] <- PercentageFeatureSet(x.seurat, pattern = "^HB[^(P)]")
  hb.genes <- rownames(x.seurat)[grep("^HB[^(P)]",rownames(x.seurat))]
  
  ###### FUNCTIONS BLOCK ########
  # QC plotting function
  QC_plotting <- function(x.QC, stage, mt.cutoff, ribo.cutoff, hb.cutoff, nFeature.min, nFeature.max) {
    
    # Visualize QC metrics as a violin plot
    VlnPlot(x.QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,".pdf"))
    
    VlnPlot(x.QC, features = "percent.mt", pt.size = 0.1) + 
      geom_hline(aes(yintercept = mt.cutoff), linetype = "dashed")
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_percent-mt.pdf"))
    
    VlnPlot(x.QC, features = "percent.ribo", pt.size = 0.1) + 
      geom_hline(aes(yintercept = ribo.cutoff), linetype = "dashed")
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_percent_ribo.pdf"))
    
    VlnPlot(x.QC, features = "percent.hb", pt.size = 0.1) + 
      geom_hline(aes(yintercept = hb.cutoff), linetype = "dashed")
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_percent_hb.pdf"))
    
    VlnPlot(x.QC, features = "nFeature_RNA", pt.size = 0.1) + 
      geom_hline(aes(yintercept = nFeature.min), linetype = "dashed")
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_nFeatures.pdf"))
    
    VlnPlot(x.QC, features = "nCount_RNA", pt.size = 0.1)
    ggsave2(paste0(Output.dir, "Vlnplot_QC_", stage,"_nCount_RNA.pdf"))
    
    # FeatureScatter plot nCount_RNA vs. percent.mt & nFeature_RNA
    plot1 <- FeatureScatter(x.QC, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(x.QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    CombinePlots(plots = list(plot1, plot2))
    ggsave2(paste0(Output.dir, "Scatterplot_QC_", stage,".pdf"))
    
    # Make one large figure of all QC-features
    feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb")
    VlnPlot(x.QC, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) +
      NoLegend()
    ggsave2(paste0(Output.dir, "All_QC_feats_", stage,".pdf"))
  }
  
  # Function for scTransform normalization
  scTransform_Normalization <- function(x,  regress_vars = FALSE, stage = "") {
    
    if (regress_vars == TRUE) {
      x = SCTransform(x, vars.to.regress = c("nFeature_RNA", "percent.mt"), verbose = T)
    } else if (regress_vars == FALSE) {
      x = SCTransform(x, verbose = T)
    }
    
    # plot variable features with and without labels
    top10_var <- head(VariableFeatures(x), 10)
    
    plot1 = VariableFeaturePlot(x)
    ggsave2(plot = plot1, filename = paste0(Output.dir, stage, "_scTransform_var_plot.pdf"))
    plot2 = LabelPoints(plot = plot1, points = top10_var, repel = TRUE)
    ggsave2(plot = plot2, filename = paste0(Output.dir, stage, "_scTransform_var_plot_top_10.pdf"))
    #var_plot = plot1 + plot2
    #ggsave2(plot = var_plot, filename = paste0(Output.dir, stage, "_scTransform_var_plot_combined.pdf"))
    
    # Do PCA on scTranformed object
    x = RunPCA(x)
    
    # Plotting PCA results
    PCA_plot = DimPlot(x, reduction = "pca")
    ggsave2(plot = PCA_plot, filename = paste0(Output.dir, stage, "_scTransform_PCA_plot.pdf"))
    
    PCA_heatmap = DimHeatmap(x, dims = 1:15, cells = 500, balanced = TRUE)
    ggsave2(plot = PCA_heatmap, filename = paste0(Output.dir, stage, "_scTransform_PCA_heatmap.pdf"))
    
    # Determine he dimensionality of the data
    elbow_plot = ElbowPlot(x, reduction = "pca", ndims = 25)
    ggsave2(plot = elbow_plot, filename = paste0(Output.dir, stage, "_scTransform_elbow_plot.pdf"))
    
    # Do UMAP on seurat object
    x = RunUMAP(x, dims = 1:30)
    
    # Run clustering on the object
    x = FindNeighbors(object = x, reduction = "pca", dims = 1:30, verbose = FALSE)
    x = FindClusters(object = x, resolution = 0.7, verbose = FALSE)
    
    # Plotting UMAP
    UMAP_plot = DimPlot(x, label = TRUE, repel = TRUE)
    ggsave2(plot = UMAP_plot, filename = paste0(Output.dir, stage, "_scTransform_UMAP.pdf"))
    
    return(x)
  }
  
  # Function for Log2 normalisation
  Log2_Normalization <- function(x, regress_vars = FALSE, stage = "") {
    x = NormalizeData(x, normalization.method = "LogNormalize")
    
    x = FindVariableFeatures(x, selection.method = "vst", verbose = T)
    top10_var <- head(VariableFeatures(x), 10)
    
    # plot variable features with and without labels
    plot1 = VariableFeaturePlot(x)
    ggsave2(plot = plot1, filename = paste0(Output.dir, stage, "_Log2_var_plot.pdf"))
    plot2 = LabelPoints(plot = plot1, points = top10_var, repel = TRUE)
    ggsave2(plot = plot2, filename = paste0(Output.dir, stage, "_Log2_var_plot_top_10.pdf"))
    var_plot = plot1 + plot2
    ggsave2(plot = var_plot, filename = paste0(Output.dir, stage, "_Log2_var_plot_combined.pdf"))
    
    all.genes = rownames(x)
    if (regress_vars == TRUE) {
      x = ScaleData(x, vars.to.regress = c("nFeature_RNA", "percent.mt"),
                    verbose = T, features = all.genes)
    } else if (regress_vars == FALSE) {
      x = ScaleData(x, verbose = T, features = all.genes)
    }
    
    # Do PCA on nromalised data data
    x = RunPCA(x, features = VariableFeatures(object = x), verbose = T)
    
    # Plotting PCA results
    PCA_plot = DimPlot(x, reduction = "pca")
    ggsave2(plot = PCA_plot, filename = paste0(Output.dir, stage, "_Log2_PCA_plot.pdf"))
    
    PCA_heatmap = DimHeatmap(x, dims = 1:15, cells = 500, balanced = TRUE)
    ggsave2(plot = PCA_heatmap, filename = paste0(Output.dir, stage, "_Log2_PCA_heatmap.pdf"))
    
    # Determine he dimensionality of the data
    elbow_plot = ElbowPlot(x, reduction = "pca", ndims = 25)
    ggsave2(plot = elbow_plot, filename = paste0(Output.dir, stage, "_Log2_elbow_plot.pdf"))
    
    x = RunUMAP(x, dims = 1:15, verbose = T)
    
    # Plotting UMAP
    UMAP_plot = DimPlot(x)
    ggsave2(plot = UMAP_plot, filename = paste0(Output.dir, stage, "_Log2_UMAP.pdf"))
    return(x)
  }
  ##### END OF FUNCTION BLOCK #########
  
  # Plotting prefiltering
  Plotting_QC = QC_plotting(x.QC = x.seurat, stage = "prefiltering", mt.cutoff = mt.cutoff,
                            ribo.cutoff = ribo.cutoff, hb.cutoff = hb.cutoff,
                            nFeature.min = nFeature.min, nFeature.max = nFeature.max)
  
  #### Before running scDblFinder for doublet/multiplet removal, do initial QC filtering:
  
  # Subset cells with low genes and genes that are expressed in at least 3 cells
  selected_cells = WhichCells(x.seurat, expression = nFeature_RNA > nFeature.min)
  selected_features = rownames(x.seurat) [Matrix::rowSums(x.seurat) > 3]
  
  x.seurat.filter.min = subset(x.seurat, features = selected_features, cells = selected_cells)
  dim(x.seurat)
  dim(x.seurat.filter.min)
  
  # Plot object after initial filtering
  Plotting_QC = QC_plotting(x.QC = x.seurat.filter.min, stage = "postfiltering_empty_droplets", mt.cutoff = mt.cutoff,
                            ribo.cutoff = ribo.cutoff, hb.cutoff = hb.cutoff,
                            nFeature.min = nFeature.min, nFeature.max = nFeature.max)
  
  ### Run scDblFinder to remove any doublets or multiplets from the dataset
  # Normalise the data either with Log2 or scTransform
  if (log2_norm_flag == TRUE & scT_norm_flag == FALSE) {
    print("Normalizing with log2")
    norm.seurat.filter.min = Log2_Normalization(x = x.seurat.filter.min, regress_vars = FALSE, stage = "postfiltering_empty_droplets_Log2")
  } else if (scT_norm_flag == TRUE & log2_norm_flag == FALSE) {
    print("Normalizing with SCTransform")
    norm.seurat.filter.min = scTransform_Normalization(x = x.seurat.filter.min, regress_vars = FALSE, stage = "postfiltering_empty_droplets_scT")
  }
  
  # Convert to sce object and run scDbFinder before reconverting
  doublet.sce = as.SingleCellExperiment(norm.seurat.filter.min)
  doublet.sce = scDblFinder(doublet.sce)
  norm.seurat.filter.min = as.Seurat(doublet.sce)
  
  # Visualise identified singlets/doublets with a UMAP
  scDblFinder_UMAP = DimPlot(norm.seurat.filter.min, group.by = ("scDblFinder.class"))
  ggsave2(plot = scDblFinder_UMAP, filename = paste0(Output.dir, "ScDblFinder_UMAP.pdf"))
  
  scDblFinder_VlnPlot = VlnPlot(norm.seurat.filter.min, features = "nFeature_RNA", group.by = "scDblFinder.class", pt.size = 0.1)
  ggsave2(plot = scDblFinder_VlnPlot, filename = paste0(Output.dir, "ScDblFinder_VlnPlot_nFeature.pdf"))
  
  scDblFinder_VlnPlot = VlnPlot(norm.seurat.filter.min, features = "nCount_RNA", group.by = "scDblFinder.class", pt.size = 0.1)
  ggsave2(plot = scDblFinder_VlnPlot, filename = paste0(Output.dir, "ScDblFinder_VlnPlot_nCount.pdf"))
  
  scDblFinder_VlnPlot = VlnPlot(norm.seurat.filter.min, features = "percent.mt", group.by = "scDblFinder.class", pt.size = 0.1)
  ggsave2(plot = scDblFinder_VlnPlot, filename = paste0(Output.dir, "ScDblFinder_VlnPlot_percent_mt.pdf"))
  
  scDblFinder_VlnPlot = VlnPlot(norm.seurat.filter.min, features = "percent.ribo", group.by = "scDblFinder.class", pt.size = 0.1)
  ggsave2(plot = scDblFinder_VlnPlot, filename = paste0(Output.dir, "ScDblFinder_VlnPlot_percent_ribo.pdf"))
  
  scDblFinder_VlnPlot = VlnPlot(norm.seurat.filter.min, features = "percent.hb", group.by = "scDblFinder.class", pt.size = 0.1)
  ggsave2(plot = scDblFinder_VlnPlot, filename = paste0(Output.dir, "ScDblFinder_VlnPlot_percent_hb.pdf"))
  
  # Filter away doublets from original seurat object to only singlets
  x.seurat.filter.min$scDblFinder.class = norm.seurat.filter.min$scDblFinder.class
  x.seurat.filter.singlets <- subset(x.seurat.filter.min, subset = scDblFinder.class == "singlet")
  
  # Plot object after scDblFinder
  Plotting_QC = QC_plotting(x.QC = x.seurat.filter.singlets, stage = "post_scDblFinder", mt.cutoff = mt.cutoff,
                            ribo.cutoff = ribo.cutoff, hb.cutoff = hb.cutoff,
                            nFeature.min = nFeature.min, nFeature.max = nFeature.max)
  
  
  # Filtering after scDblFinder doublet finder
  selected_mt.ribo.hb = WhichCells(x.seurat.filter.singlets, expression = percent.mt < mt.cutoff &
                                     percent.ribo < ribo.cutoff &
                                     percent.hb < hb.cutoff)
  
  x.seurat.filter.max = subset(x.seurat.filter.singlets, cells = selected_mt.ribo.hb)
  dim(x.seurat)
  dim(x.seurat.filter.min)
  dim(x.seurat.filter.singlets)
  dim(x.seurat.filter.max)
  
  # Remove mitochondrial genes
  x.seurat.filter.max <- x.seurat.filter.max[!grepl("^MT-", rownames(x.seurat.filter.max)), ]	

  # Plot object after scDblFinder and filtering mt, ribo and hb
  Plotting_QC = QC_plotting(x.QC = x.seurat.filter.max, stage = "post_scDblFinder_&_filtering", mt.cutoff = mt.cutoff,
                            ribo.cutoff = ribo.cutoff, hb.cutoff = hb.cutoff,
                            nFeature.min = nFeature.min, nFeature.max = nFeature.max)
  
  # Normalise the data either with Log2 or scTransform
  if (log2_norm_flag == TRUE & scT_norm_flag == FALSE) {
    norm.seurat.filter.max = Log2_Normalization(x = x.seurat.filter.max, regress_vars = FALSE, stage = "post_scDblFinder_&_filtering_Log2")
  } else if (scT_norm_flag == TRUE & log2_norm_flag == FALSE) {
    norm.seurat.filter.max = scTransform_Normalization(x = x.seurat.filter.max, regress_vars = FALSE, stage = "post_scDblFinder_&_filtering_scT")
  }
  
  #### Calculate cell-cycle scores ####
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  norm.seurat.filter.max = CellCycleScoring(norm.seurat.filter.max, s.features = s.genes, 
                                            g2m.features = g2m.genes, set.ident = TRUE)
  
  # Plot cell cycle phases
  Phase_UMAP = DimPlot(norm.seurat.filter.max, group.by = "Phase")
  ggsave2(plot = Phase_UMAP, filename = paste0(Output.dir, "Cell-Cycle_Phase_UMAP.pdf"))
  
  # Visualize the distribution of cell cycle markers across
  VlnPlot(norm.seurat.filter.max, features = c("S.Score","G2M.Score"), group.by = "seurat_clusters", pt.size = 0.1)
  ggsave2(paste0(Output.dir, "Cell.Cycle_Phase_VlnPlot.pdf"))
  
  ## Save filtered and clustered data as new seurat object ##
  SaveH5Seurat(norm.seurat.filter.max, file = paste0(Output.dir, x, "_Filtered_Endometrium.h5seurat"), overwrite = TRUE)
  
  # Plotting additional UMAPs
  FeaturePlot(norm.seurat.filter.max, features = "percent.mt")
  ggsave2(paste0(Output.dir, "UMAP_percent_mt.pdf"))
  
  FeaturePlot(norm.seurat.filter.max, features = "percent.ribo")
  ggsave2(paste0(Output.dir, "UMAP_percent_ribo.pdf"))
  
  FeaturePlot(norm.seurat.filter.max, features = "percent.hb")
  ggsave2(paste0(Output.dir, "UMAP_percent_hb.pdf"))
  
  FeaturePlot(norm.seurat.filter.max, features = "nFeature_SCT")
  ggsave2(paste0(Output.dir, "UMAP_nFeature_SCT.pdf"))
  
  FeaturePlot(norm.seurat.filter.max, features = "nFeature_RNA")
  ggsave2(paste0(Output.dir, "UMAP_nFeature_RNA.pdf"))
  
  FeaturePlot(norm.seurat.filter.max, features = "nCount_SCT")
  ggsave2(paste0(Output.dir, "UMAP_nCount_SCT.pdf"))
  
  FeaturePlot(norm.seurat.filter.max, features = "nCount_RNA")
  ggsave2(paste0(Output.dir, "UMAP_nCount_RNA.pdf"))
  
  FeaturePlot(norm.seurat.filter.max, features = "S.Score")
  ggsave2(paste0(Output.dir, "UMAP_S_score.pdf"))
  
  FeaturePlot(norm.seurat.filter.max, features = "G2M.Score")
  ggsave2(paste0(Output.dir, "UMAP_G2M_Score.pdf"))
  
}

# Apply pre-processing function on all samples
Output.list = lapply(Input.list, pre_processing.fun)
