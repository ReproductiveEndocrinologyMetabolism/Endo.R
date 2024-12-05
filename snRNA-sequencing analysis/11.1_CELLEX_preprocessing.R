#!/usr/bin/env Rscript
library(dplyr)

# Function to preprocess the Seurat object to a gene count matrix and metadata for CELLECT analysis
CELLEX_preprocessing <- function(Project.name, seurat.path, seurat.ident.main, seurat.ident.sub,
                                 generate.matrix.step = TRUE, ensembl.conversion = TRUE, format.matrix.step = FALSE) {
  
  # Generate or set the output directory
  if (dir.exists(path = "Output/CELLEX_preprocessing") == FALSE) {
    print("Output/CELLEX_preprocessing")
    dir.create(path = "Output/CELLEX_preprocessing", recursive = TRUE)
    Output.dir = "Output/CELLEX_preprocessing/"
  } else if (dir.exists(path = "Output/CELLEX_preprocessing") == TRUE) {
    print("Directory exists")
    Output.dir = "Output/CELLEX_preprocessing/"
  } else {
    print("Error with output directory")
  }
  
  if (generate.matrix.step == TRUE && format.matrix.step == FALSE) {
    
    print("Generating the matrices")
    
    # Load the packages
    library(Seurat)
    library(SeuratDisk)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    
    # Load the Seurat object
    seurat.x = LoadH5Seurat(seurat.path, assays = "RNA")
    
    # Generate a dataframe with the count data
    genecount.df = data.frame(GetAssayData(object = seurat.x, slot = "counts", assay = "RNA"))
    
    if (ensembl.conversion == TRUE) {
      
      print("Saving the matrix with symbols")
      
      # Save the genecount matrix with symbols to .csv
      write.csv(x = genecount.df, file = paste0(Output.dir, Project.name, "_genecount_data_Symbol.csv"), row.names = TRUE)
      test.df = read.csv(file = paste0(Output.dir, Project.name, "_genecount_data_Symbol.csv"))
      print("Check the SYMBOL .csv")
      head(test.df)
      test.df = NULL
      
      print("Converting to Ensembl")
      
      # Convert symbol to ensembl ID
      genecount.df$cell_id = mapIds(org.Hs.eg.db,
                                    keys=rownames(genecount.df), 
                                    column="ENSEMBL",
                                    keytype="SYMBOL",
                                    multiVals="first")
      
      # Add the new Ensembl id column first
      genecount.df <- genecount.df[, c(ncol(genecount.df), 1:(ncol(genecount.df)-1))]
      
      print("Removing NA's")
      
      # Remove rows with Ensembl ID that are NA
      genecount.df = subset(genecount.df, !is.na(cell_id))
      
      # Rename the cell_id column to Ensembl
      colnames(genecount.df)[colnames(genecount.df) == "cell_id"] <- "Ensembl_id"
      
      print("Saving the matrix with Ensembl")
      
      # Save the genecount matrix to .csv
      write.csv(x = genecount.df, file = paste0(Output.dir, Project.name, "_genecount_data_Ensembl.csv"), row.names =FALSE)
      test.df = read.csv(file = paste0(Output.dir, Project.name, "_genecount_data_Ensembl.csv"))
      
      print("Check the ENSEMBL .csv")
      head(test.df)
      test.df = NULL
      
      genecount.df = NULL
    } else if (ensembl.conversion == FALSE) {
      # Save the genecount matrix to .csv
      write.csv(x = genecount.df, file = paste0(Output.dir, Project.name, "_genecount_data_Symbol.csv"), row.names =FALSE)
    }
    
    # Remove the dataframe to reduce memory load
    genecount.df = NULL
    
    # Generate the metadata with the cell subtypes and save.
    data.x = data.frame(cell_id = rownames(seurat.x@meta.data), cell_type = seurat.x[[seurat.ident.main]], 
                        group = seurat.x$Group_Stage, 
                        sample = seurat.x$orig.ident)
    
    print("Metadata matrix of main celltypes")
    print(head(data.x))
    write.csv(data.x, file = paste0(Output.dir, Project.name, "_metadata_Main.csv"), row.names = FALSE)
    
    # Generate the metadata with the cell subtypes and save.
    data.x = data.frame(cell_id = rownames(seurat.x@meta.data), cell_type = seurat.x[[seurat.ident.sub]], 
                        group = seurat.x$Group_Stage, 
                        sample = seurat.x$orig.ident)
    
    print("Metadata matrix of subtypes")
    print(head(data.x))
    write.csv(data.x, file = paste0(Output.dir, Project.name, "_metadata_Sub.csv"), row.names = FALSE)
    
    # Generate the metadata with the cell subtypes and save.
    data.x = data.frame(cell_id = rownames(seurat.x@meta.data), 
                        main_celltype = seurat.x[[seurat.ident.main]],
                        sub_celltype = seurat.x[[seurat.ident.sub]],
                        group = seurat.x$Group_Stage, 
                        sample = seurat.x$orig.ident)
    
    print("Metadata matrix of subtypes")
    print(head(data.x))
    write.csv(data.x, file = paste0(Output.dir, Project.name, "_metadata.csv"), row.names = FALSE)
    
    print("Done with generating matrices")
    
  } else if (generate.matrix.step == FALSE && format.matrix.step == TRUE) {
    
    print("Formatting the matrices")
    
    # Load required packages
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    
    print("Loading the .csv file")
    
    # Load the genecount matrix
    gene.csv = read.csv(file = paste0(Output.dir, Project.name, "_genecount_data_Ensembl.csv"))
    
    # Convert symbol to ensembl
    gene.csv$cell_id = mapIds(org.Hs.eg.db,
                              keys=gene.csv$cell_symbol, 
                              column="ENSEMBL",
                              keytype="SYMBOL",
                              multiVals="first")
    
    # Check how many that are not mapped
    print(paste("Unmapped genes:", sum(is.na(gene.csv$cell_id))))
    gene.na = gene.csv[is.na(gene.csv$cell_id),]
    head(gene.na$cell_symbol)
    
    # Create the new genecount with Ensembl
    gene.csv <- gene.csv[, c(ncol(gene.csv), 1:(ncol(gene.csv)-1))]
    
    print("Saving the merged symbol and ensembl matrix")
    
    # Save csv with both symbol and Ensembl
    write.csv(x = gene.csv, file = paste0(Output.dir, Project.name, "_genecount_data_Symbol_Ensembl.csv"), row.names =FALSE)
    
    # Remove the cell_symbol column
    gene.csv$cell_symbol = NULL
    
    # Remove rows with Ensembl ID that are NA
    gene.csv = subset(gene.csv, !is.na(cell_id))
    
    print("Saving the ensembl matrix")
    
    # Save csv with only Ensembl ID
    write.csv(x = gene.csv, file = paste0(Output.dir, Project.name, "_genecount_data_Ensembl.csv"), row.names =FALSE)
    
  }
  
}

# Set to whether run CELLEX preprocessing on the group data or celltype
do.group.CELLECT = TRUE
do.celltype.CELLECT = FALSE

# Run the CELLEX preprocessing function on either group or major celltypes
if (do.group.CELLECT == TRUE) {
  CELLEX_Control = CELLEX_preprocessing(Project.name = "Endo_Combined_Control",
                                        seurat.path = "Data/Seurat_Subsets/Endo_Combined_Control_subset.h5seurat",
                                        seurat.ident.main = "Labelled_Clusters_SCT.1",
                                        seurat.ident.sub = "Combined_labels")
  CELLEX_PCOS = CELLEX_preprocessing(Project.name = "Endo_Combined_PCOS",
                                        seurat.path = "Data/Seurat_Subsets/Endo_Combined_PCOS_W0_subset.h5seurat",
                                        seurat.ident.main = "Labelled_Clusters_SCT.1",
                                        seurat.ident.sub = "Combined_labels")
  CELLEX_Baseline = CELLEX_preprocessing(Project.name = "Endo_Combined_Baseline",
                                        seurat.path = "Data/Seurat_Subsets/Endo_Combined_Control_PCOS_W0_subset.h5seurat",
                                        seurat.ident.main = "Labelled_Clusters_SCT.1",
                                        seurat.ident.sub = "Combined_labels")
} else if (do.celltype.CELLECT == TRUE) {
  CELLEX_epithelium = CELLEX_preprocessing(Project.name = "Endo_10x_Epithelium",
                                           seurat.path = "Output/5_Epithelium_Clustering/Endo_All_Epithelium_reclustered_labelled.h5seurat",
                                           seurat.ident.main = "Labelled_Clusters_SCT.1",
                                           seurat.ident.sub = "Epithelium_labelled")
  
  CELLEX_Stroma = CELLEX_preprocessing(Project.name = "Endo_10x_Stroma",
                                       seurat.path = "Output/5_Stromal_Clustering/Endo_All_Stromal_reclustered_labelled.h5seurat",
                                       seurat.ident.main = "Labelled_Clusters_SCT.1",
                                       seurat.ident.sub = "Stromal_labelled")
  
  CELLEX_immune = CELLEX_preprocessing(Project.name = "Endo_10x_Immune",
                                       seurat.path = "Output/5_Immune_Clustering/Endo_All_Immune_reclustered_labelled.h5seurat",
                                       seurat.ident.main = "Labelled_Clusters_SCT.1",
                                       seurat.ident.sub = "Immune_labelled")
  
  CELLEX_endothelium = CELLEX_preprocessing(Project.name = "Endo_10x_Endothelial",
                                            seurat.path = "Output/5_Endothelial_Clustering/Endo_All_Endothelial_reclustered_labelled.h5seurat",
                                            seurat.ident.main = "Labelled_Clusters_SCT.1",
                                            seurat.ident.sub = "Endothelial_labelled")
}
