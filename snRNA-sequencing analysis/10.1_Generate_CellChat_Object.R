#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(CellChat)
library(openxlsx)
library(patchwork)

# Generate or set the output directory
if (dir.exists(path = "Output/9_CellChat") == FALSE) {
  print("Output/9_CellChat")
  dir.create(path = "Output/9_CellChat", recursive = TRUE)
  Output.dir = "Output/9_CellChat/"
} else if (dir.exists(path = "Output/9_CellChat") == TRUE) {
  print("Directory exists")
  Output.dir = "Output/9_CellChat/"
} else {
  print("Error with output directory")
}

Project_name = "Endo_All"
select.idents = "Combined_labels"
group.idents = "Group_Stage"

# Add idents here to be kept in in select.idents. Used to remove Undefined celltype
levels.idents = c("Lumenal", "SOX9+ LGR5+", "SOX9+ LGR5-",
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
Norm.assay = "RNA"
load.seurat.subset = TRUE
normalize.seurat = TRUE
group.pattern = ".*Endo_All_(.*?)_Combined_labels.h5seurat"

if (load.seurat.subset == TRUE) {
  
  # List group subsetted seurats to be loaded
  seurat.groups = list.files("Output/9_Main_relabelled", 
                             pattern = "Endo_All_.*_Combined_labels.h5seurat",
                             full.names = TRUE)
  
} else if (load.seurat.subset == FALSE) {
  
  # Loading the seurat object
  seurat.x = LoadH5Seurat(file = "Output/9_Main_relabelled/Endo_All_Combined_labels.h5seurat")
  
  # Setting up the Seurat object
  Idents(object = seurat.x) =  group.idents
  DefaultAssay(seurat.x) = Norm.assay
  seurat.groups = levels(seurat.x)
  
  if (normalize.seurat == TRUE) {
    seurat.subset = NormalizeData(seurat.subset)
  }
  
}

CellChat.subset <- function(seurat.subset) {
  
  # Generating the CellChat object
  CC.input <- GetAssayData(seurat.subset, assay = "RNA", slot = "data") # Normalized data matrix
  #Idents(object = seurat.subset) = select.idents
  CC.labels <- Idents(seurat.subset)
  CC.meta <- data.frame(labels = CC.labels, row.names = names(CC.labels)) # Create a dataframe of the cell labels
  CC.x <- createCellChat(object = CC.input, meta = CC.meta, group.by = "labels")
  
  # Set the CellChat database
  CC.x@DB <- CellChatDB.human
  
  # Preprocessing the CellChat object before cell-cell communication analysis
  CC.x <- subsetData(object = CC.x) # Perhaps only subset DEG genes in the future
  CC.x <- identifyOverExpressedGenes(CC.x)
  CC.x <- identifyOverExpressedInteractions(CC.x)
  CC.x <- projectData(CC.x, PPI.human) # Projecting the data on protein-protein interaction (PPI) network
  
  # Inference of cell-cell communication network
  CC.x <- computeCommunProb(CC.x, population.size = TRUE)
  CC.x <- filterCommunication(CC.x, min.cells = 10)
  
  # Inferring the cell-cell communication pathway signaling
  CC.x <- computeCommunProbPathway(CC.x)
  
  # Calculating the aggregated cell-cell network communication
  CC.x <- aggregateNet(CC.x)
  
  # Return the CellChat object
  return(CC.x)
  
}

#x.subset = seurat.groups[1] # FOR TESTING
CC.list = list()
for (x.subset in seurat.groups) {
  
  if (load.seurat.subset == TRUE) {
    
    # Loading the group ubsetted seurat object
    seurat.subset = LoadH5Seurat(file = x.subset)
    x.subset = gsub(group.pattern, "\\1", x.subset)
    
    if (normalize.seurat == TRUE) {
      seurat.subset = NormalizeData(seurat.subset)
    }
    
  } else if (load.seurat.subset == FALSE) {
    
    # Subset the seurat object and extract the group
    seurat.subset = subset(seurat.x, idents = x.subset)
    
  }
  
  print(paste("CellChat is running on", x.subset))
  
  # Setting up the Seurat object
  DefaultAssay(seurat.subset) = Norm.assay
  Idents(object = seurat.subset) = select.idents
  
  if (is.null(levels.idents) == FALSE) {
    print("Setting the idents")
    seurat.subset = subset(seurat.subset, idents = levels.idents)
    seurat.subset@active.ident = factor(x = seurat.subset@active.ident, levels = levels.idents)
  }

  # Run the CellChat subsetting function
  CC.save = CellChat.subset(seurat.subset = seurat.subset)
  
  # Save the CellChat object
  saveRDS(CC.save, file = paste0(Output.dir, Project_name, "_CellChat_object_", x.subset, ".rds"))
  
  # Extract the data to a dataframe and save it as an excel
  
  # Append the CC.save to a list
  CC.list[x.subset] = CC.save
 
}

#CC.list = list()
#CC.list["Control"] = readRDS(file = "Output/9_CellChat/Endo_All_CellChat_object_Control.rds")
#CC.list["PCOS_W0"] = readRDS(file = "Output/9_CellChat/Endo_All_CellChat_object_PCOS_W0.rds")
#CC.list["PCOS_W16_LS"] = readRDS(file = "Output/9_CellChat/Endo_All_CellChat_object_PCOS_W16_LS.rds")
#CC.list["PCOS_W16_Met"] = readRDS(file = "Output/9_CellChat/Endo_All_CellChat_object_PCOS_W16_Met.rds")

# Merge the CellChat objects
CC.baseline = mergeCellChat(CC.list[c(1,2)], add.names = c("Control", "PCOS_W0"))
CC.metformin = mergeCellChat(CC.list[c(2,4)], add.names = c("PCOS_W0", "PCOS_W16_Met"))
CC.lifestyle = mergeCellChat(CC.list[c(2,3)], add.names = c("PCOS_W0", "PCOS_W16_LS"))
CC.combined <- mergeCellChat(CC.list[c(1, 2, 4, 3)], add.names = names(CC.list)[c(1, 2, 4, 3)])

# Save the merged CellChat object
saveRDS(CC.baseline, file = paste0(Output.dir, Project_name, "_CellChat_Merged_object_Control_PCOS.rds"))
saveRDS(CC.metformin, file = paste0(Output.dir, Project_name, "_CellChat_Merged_object_PCOS_Metformin.rds"))
saveRDS(CC.lifestyle, file = paste0(Output.dir, Project_name, "_CellChat_Merged_object_PCOS_Lifestyle.rds"))
saveRDS(CC.combined, file = paste0(Output.dir, Project_name, "_CellChat_Merged_object_All.rds"))
