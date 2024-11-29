#!/usr/bin/env Rscript

print("Loading packages")
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(future)

### Flags and directories to be set ###
LogNorm.flag = FALSE
SCT.flag = TRUE
sample.reference.flag = FALSE
do.integration.flag = TRUE
load.anchors = TRUE
integration.feat = 2500

Project_name = "Endo_All"
sample.integration = list.files("Output/0_Pre-processed/")
sample.references = c(1:5) # These are the controls
Input.dir = "Output/0_Pre-processed/"
anchor.dir = "Output/1_Integrated/Endo_All_endo_anchors.RDS"

# Creating the output directory
if (dir.exists(path = paste0("Output/1_Integrated/")) == FALSE) {
  print(paste0("Generating output directory ", "1_Integrated"))
  dir.create(path = paste0("Output/1_Integrated/"), recursive = TRUE)
  Output.dir = paste0("Output/1_Integrated/")
} else if (dir.exists(path = paste0("Output/1_Integrated/")) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/1_Integrated/")
} else {
  print("Error with output directory")
}

### Generating the Seurat object ###

# Load integrated seurat object if it exists
if (length(list.files(Output.dir, pattern = "integrated.h5seurat", full.names = TRUE)) == 1) {
  print(paste0("Loading ", list.files(Output.dir, pattern = "_integrated.h5seurat", full.names = FALSE)))
  endo.integrated = LoadH5Seurat(file = list.files(Output.dir, pattern = "_integrated.h5seurat", full.names = TRUE))
} else if (length(list.files(Output.dir, pattern = "integrated.h5seurat", full.names = TRUE)) == 0 & load.anchors == FALSE) {
  # Perform integration if the integration has not been performed before
  print("Generating integrated seurat object")
  
  # Loop to load sample seurat objects
  seurat.dir = c()
  for (x in sample.integration) {
    sample.dir = x
    seurat.path = paste0(Input.dir, sample.dir, "/", sample.dir, "_Filtered_Endometrium.h5seurat")
    seurat.dir = c(seurat.dir, seurat.path)
  }
  seurat.list = lapply(seurat.dir, function(x) LoadH5Seurat(file = x))
  
} else if (load.anchors == TRUE) {
  print("Anchors are to be loaded!")
}

if (LogNorm.flag == TRUE & load.anchors == FALSE) {
  print("Log Normalised")
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = TRUE)
    x <- FindVariableFeatures(x, verbose = TRUE)
  })
  # Finding integration anchors and integrating
  if (sample.reference.flag == TRUE) {
    print("Using reference")
    endo.anchors <- FindIntegrationAnchors(object.list = seurat.list,
                                           normalization.method = "LogNormalize",
                                           reference = sample.references, 
                                           reduction = "rpca", dims = 1:50)
    
    # The anchor dataset is saved so that it can be used on the cluster
    saveRDS(endo.anchors, paste0(Output.dir,Project_name, "_endo_anchors.RDS"))
    
    
  } else if (sample.reference.flag == FALSE) {
    print("Using reference")
    endo.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                           reduction = "rpca", dims = 1:50)
    endo.integrated <- IntegrateData(anchorset = endo.anchors, dims = 1:50)
    endo.integrated <- ScaleData(endo.integrated, verbose = FALSE)
    endo.integrated <- RunPCA(endo.integrated, verbose = FALSE)
    endo.integrated <- RunUMAP(endo.integrated, dims = 1:50)
  }
  
  # The anchor dataset is saved
  saveRDS(endo.anchors, paste0(Output.dir,Project_name, "_endo_anchors.RDS"))
  
  # The data is integrated and saved
  if (do.integration.flag == TRUE) {
    endo.integrated <- IntegrateData(anchorset = endo.anchors, normalization.method = "LogNormalize", dims = 1:50)
    endo.integrated <- ScaleData(endo.integrated, verbose = FALSE)
    endo.integrated <- RunPCA(endo.integrated, verbose = FALSE)
    endo.integrated <- RunUMAP(endo.integrated, dims = 1:50)
    SaveH5Seurat(endo.integrated, paste0(Output.dir,Project_name, "_integrated.h5seurat"), overwrite = TRUE)
  }
  
} else if (SCT.flag == TRUE & load.anchors == FALSE) {
  print("Integration of SCT transformed object")
  # Preparing SC transformed objects for integration
  #seurat.list <- lapply(X = seurat.list, FUN = SCTransform) ### They have already been normalised SCTransformed
  features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = integration.feat)
  seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
  seurat.list <- lapply(X = seurat.list, FUN = RunPCA, features = features)
  
  # Finding integration anchors and integrating
  if (sample.reference.flag == TRUE) {
    print("Using reference")
    endo.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                           normalization.method = "SCT", 
                                           reference = sample.references, 
                                           anchor.features = features,
                                           dims = 1:30,
                                           reduction = "rpca", 
                                           k.anchor = 5)
  } else if (sample.reference.flag == FALSE) {
    endo.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                           normalization.method = "SCT", 
                                           anchor.features = features,
                                           dims = 1:30,
                                           reduction = "rpca", 
                                           k.anchor = 5)
  }
  
  # The anchor dataset is saved
  saveRDS(endo.anchors, paste0(Output.dir,Project_name, "_endo_anchors.RDS"))
  
  # The data is integrated and saved
  if (do.integration.flag == TRUE) {
    endo.integrated <- IntegrateData(anchorset = endo.anchors, normalization.method = "SCT", dims = 1:30)
    endo.integrated <- RunPCA(endo.integrated, verbose = TRUE)
    endo.integrated <- RunUMAP(endo.integrated, reduction = "pca", dims = 1:30)
    endo.integrated <- FindNeighbors(endo.integrated, dims = 1:30)
    endo.integrated <- FindClusters(endo.integrated)
  }
  
} else if (load.anchors == TRUE & do.integration.flag == TRUE) {
  print("Loading the anchors!")
  endo.anchors = readRDS(file = anchor.dir)
  
  # The data is integrated and saved
  if (do.integration.flag == TRUE) {
    
    if (SCT.flag == TRUE) {
      print("Integration with SCT and loaded anchors")
      endo.integrated <- IntegrateData(anchorset = endo.anchors, normalization.method = "SCT", dims = 1:30)
      endo.integrated <- RunPCA(endo.integrated, verbose = TRUE)
      endo.integrated <- RunUMAP(endo.integrated, reduction = "pca", dims = 1:30)
      endo.integrated <- FindNeighbors(endo.integrated, dims = 1:30)
      endo.integrated <- FindClusters(endo.integrated)
    } else if (LogNorm.flag == TRUE) {
      print("Integration with log normalisation and loaded anchors")
      endo.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                             reduction = "rpca", dims = 1:50)
      endo.integrated <- IntegrateData(anchorset = endo.anchors, dims = 1:50)
      endo.integrated <- ScaleData(endo.integrated, verbose = FALSE)
      endo.integrated <- RunPCA(endo.integrated, verbose = FALSE)
      endo.integrated <- RunUMAP(endo.integrated, dims = 1:50)
    }
  }
  
}

# Saving the seurat object
SaveH5Seurat(object = endo.integrated, filename = paste0(Output.dir,Project_name,"_integrated.h5seurat"), overwrite = TRUE)

# Plotting the integrated object
print("Plotting of UMAP")
UMAP_Ident = DimPlot(endo.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_labelled_UMAP_Ident.pdf"), dpi = 700)
UMAP_Ident = DimPlot(endo.integrated, reduction = "umap", group.by = "orig.ident", label = TRUE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay,"_labelled_UMAP_Ident.pdf"), dpi = 700)
UMAP_Phase = DimPlot(endo.integrated, reduction = "umap", group.by = "Phase")
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay,"_UMAP_Phase.pdf"), dpi = 700)
UMAP_Label = DimPlot(endo.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay,"_labelled_UMAP_Clusters.pdf"), dpi = 700)
print("Done with basic UMAP")

# Renaming for groups and stage
endo.integrated$Group_Stage = endo.integrated$orig.ident
endo.integrated = SetIdent(endo.integrated, value = "Group_Stage")
endo.integrated = RenameIdents(endo.integrated, '1_Control' = "Control",
                               "2_Control" = "Control",
                               "3_Control" = "Control",
                               "4_Control" = "Control",
                               "5_Control" = "Control",
                               "6_PCOS_Lifestyle_W0" = "PCOS_W0",
                               "18_PCOS_Lifestyle_W16" = "PCOS_W16_LS",
                               "7_PCOS_Lifestyle_W0" = "PCOS_W0",
                               "19_PCOS_Lifestyle_W16" = "PCOS_W16_LS",
                               "8_PCOS_Lifestyle_W0" = "PCOS_W0",
                               "20_PCOS_Lifestyle_W16" = "PCOS_W16_LS",
                               "9_PCOS_EA_W0" = "PCOS_W0",
                               "10_PCOS_Metformin_W0" = "PCOS_W0",
                               "21_PCOS_Metformin_W16" = "PCOS_W16_Met",
                               "11_PCOS_Metformin_W0" = "PCOS_W0",
                               "22_PCOS_Metformin_W16" = "PCOS_W16_Met",
                               "12_PCOS_Metformin_W0" = "PCOS_W0",
                               "23_PCOS_Metformin_W16" = "PCOS_W16_Met",
                               "13_PCOS_Metformin_W0" = "PCOS_W0",
                               "24_PCOS_Metformin_W16" = "PCOS_W16_Met",
                               "14_PCOS_Metformin_W0" = "PCOS_W0",
                               "25_PCOS_Metformin_W16" = "PCOS_W16_Met",
                               "15_PCOS_Metformin_W0" = "PCOS_W0",
                               "26_PCOS_Metformin_W16" = "PCOS_W16_Met",
                               "16_PCOS_Metformin_W0" = "PCOS_W0",
                               "17_PCOS_Metformin_W0" = "PCOS_W0",
                               "27_PCOS_Metformin_W16" = "PCOS_W16_Met")
endo.integrated[["Group_Stage"]] = Idents(object = endo.integrated)
DimPlot(endo.integrated, reduction = "umap", group.by = "Group_Stage")
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay,"_UMAP_Group_Stage.pdf"), dpi = 700)

# Renaming for groups and stage
DimPlot(endo.integrated, reduction = "umap", group.by = "orig.ident")
ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay,"_UMAP_Pat-nr.pdf"), dpi = 700)

# log normalising 
DefaultAssay(endo.integrated) = "RNA"
endo.integrated = NormalizeData(endo.integrated)
endo.integrated = FindVariableFeatures(endo.integrated)
all.genes = rownames(endo.integrated)
endo.integrated = ScaleData(endo.integrated, features = all.genes)
Norm.assay = "RNA"


# Save labelled seurat object and overwrite old object
SaveH5Seurat(object = endo.integrated, filename = paste0(Output.dir,Project_name, "_", Norm.assay,"_labeled_integrated.h5seurat"), overwrite = TRUE)

# Plotting QC metrics 
Plotting_QC <- function(x, group = "") {
  
  # Generating ridgeplots
  RidgePlot(endo.integrated, features = "percent.mt", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_mtDNA.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "percent.ribo", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_ribo.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "percent.hb", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_hb.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "nCount_RNA", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_nCount_RNA.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "nFeature_RNA", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_nFeature_RNA.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "S.Score", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_S-Score.pdf"), dpi = 700)
  
  RidgePlot(endo.integrated, features = "G2M.Score", group.by = group)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Ridgeplot_G2M-Score.pdf"), dpi = 700)
  
  # Generating violin plots
  VlnPlot(endo.integrated, features = "percent.mt", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_mtDNA.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "percent.ribo", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_ribo.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "percent.hb", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_hb.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "nCount_RNA", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_nCount_RNA.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "nFeature_RNA", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_nFeature_RNA.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "S.Score", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_S-Score.pdf"), dpi = 700)
  
  VlnPlot(endo.integrated, features = "G2M.Score", group.by = group, pt.size = 0)
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_Vlnplot_G2M-Score.pdf"), dpi = 700)
  
  # Make stacked barplot of group population per cell type to determine proportionality
  # between Control and PCOS in each cell type and after treatment
  
  Idents(object = endo.integrated) <- group
  pt <- table(Idents(endo.integrated), endo.integrated$Group_Stage)
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)
  
  # Add colours for the colour blind
  barplot_prop = ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Group") +
    ylab("Proportion") +
    theme_cowplot() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, 
                                                                     vjust = 1, hjust = 1))
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_barplot_proportions_Stage.pdf"), dpi = 700)
  
  # Make stacked barplot of group population per cell type to determine proportionality
  # between all different samples
  
  Idents(object = endo.integrated) <- group
  pt <- table(Idents(endo.integrated), endo.integrated$orig.ident)
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)
  
  # Add colours for the colour blind
  barplot_prop = ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Sample") +
    ylab("Proportion") +
    theme_cowplot() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, 
                                                                     vjust = 1, hjust = 1))
  ggsave2(paste0(Output.dir,Project_name, "_", Norm.assay, "_",group, "_barplot_proportions_Sample.pdf"), dpi = 700)
  
}

Plot_CQ_Group.Stage = Plotting_QC(x = endo.integrated, group = "Group_Stage")
Plot_CQ_Pat_nr = Plotting_QC(x = endo.integrated, group = "orig.ident")
Plot_CQ_seurat_clusters = Plotting_QC(x = endo.integrated, group = "seurat_clusters")

# Make stacked barplot of group population per cell type to determine proportionality
# between Control and PCOS in each cell type and after treatment

Proportion_barplot <- function(x_group = "") {
  Idents(object = endo.integrated) <- group
  pt <- table(Idents(endo.integrated), endo.integrated[[x_group]])
  pt <- as.data.frame(pt)
  pt$Var1 <- as.character(pt$Var1)
  
  # Add colours for the colour blind
  barplot_prop = ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
    theme_bw(base_size = 15) +
    geom_col(position = "fill", width = 0.5) +
    xlab("Sample") +
    ylab("Proportion") +
    theme_cowplot() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, 
                                                                     vjust = 1, hjust = 1))
  ggsave2(plot = barplot_prop, filename = paste0(Output.dir,Project_name, "_", Norm.assay, "_", "_barplot_proportions_", x_group, ".pdf"), dpi = 700)
  
}

barplot_prop_QC = Proportion_barplot(x_group = "Group_Stage")
barplot_prop_QC = Proportion_barplot(x_group = "orig.ident")
