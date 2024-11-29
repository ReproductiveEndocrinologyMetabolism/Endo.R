#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(openxlsx)
library(stringr)
library(ggvenn)

# Depending on running on server or not, all packages cannot be loaded
server.run = FALSE #TRUE FALSE
save.subset.seurat = FALSE
skip.overlap.flag = FALSE # TRUE in endothelial

# Setting up script
#Project_name = "Endo_All_Epithelium_DGE"
celltype = "Epithelium" #"Epithelium"  #"Stromal_uSMC" # Endothelial # Immune
Project_name = paste0("Endo_All_", celltype, "_DGE")
select.idents = "Epithelium_labelled" #"Epithelium_labelled" #"Stromal_labelled"
Input.dir = paste0("Output/6_DGE_analysis_", select.idents, "/")
theme_set(theme_cowplot())
#DEG.comp = c("CtrlvsPCOS", "Overlapped", "PCOSvsMetformin", "PCOSvsLifestyle") # For GO enrichment

if (server.run == TRUE) {
  
  seurat.dir = paste0("Output/5_", celltype, "_Clustering/")
  seurat.project = paste0("Endo_All_", celltype)
  
  
} else if (server.run == FALSE) {
 
  library(enrichR)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(DOSE)
  library(pathview)
  library(msigdbr)
  
  # Extracting genelists for GO enrichment analysis
  m_Hallmark = msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, gene_symbol)

  m_Immune = msigdbr(species = "Homo sapiens", category = "C7", subcategory = "IMMUNESIGDB") %>% 
    dplyr::select(gs_name, gene_symbol)
  
  m_HPO = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "HPO") %>% 
    dplyr::select(gs_name, gene_symbol)
  
  m_msigdb_BP = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% 
    dplyr::select(gs_name, gene_symbol)

}

options(ggrepel.max.overlaps = Inf)

# Checking and/or generating output dir
if (dir.exists(path = paste0("Output/7_Selected_DEG_GO_analysis_", select.idents)) == FALSE) {
  print(paste0("Generating output directory Output/7_Selected_DEG_GO_analysis_", select.idents))
  dir.create(path = paste0("Output/7_Selected_DEG_GO_analysis_", select.idents), recursive = TRUE)
  Output.dir = paste0("Output/7_Selected_DEG_GO_analysis_", select.idents, "/")
} else if (dir.exists(path = paste0("Output/7_Selected_DEG_GO_analysis_", select.idents)) == TRUE) {
  print("Directory exists")
  Output.dir = paste0("Output/7_Selected_DEG_GO_analysis_", select.idents, "/")
} else {
  print("Error with output directory")
}

log2fc.cutoff = 0.5
q.value = 0.05
x.label = "CtrlvsPCOS"
y.label = "PCOSvsMetformin"
z.label = "PCOSvsLifestyle"
Two.cols = c("dodgerblue", "firebrick")
Three.cols = c("dodgerblue", "seashell", "firebrick")
Group.cols = c("#A0A0A0", "#D098B2", "#95BFE1", "#65D46E")
Venn.cols = c("#b89ca9", "#b3acca", "#9bb690") # Mixed colours of CtrlvsPCOS, PCOSvsMet, PCOSvsLS 

# Create DEG list of genes overlapped between DEG lists
Overlap.DEGs <- function(celltype.x) {
  
  input.dir.list = list.files(path = paste0(Input.dir, celltype.x), pattern = "*_sign_DEG.xlsx", full.names = TRUE)
  x.df.dir = input.dir.list[grep(x.label, input.dir.list)]
  y.df.dir = input.dir.list[grep(y.label, input.dir.list)]
  z.df.dir = input.dir.list[grep(z.label, input.dir.list)]
  
  df.list = list(x.df.dir, y.df.dir, z.df.dir)
  names(df.list) = c(x.label, y.label, z.label)
  
  filtered.tables = list()
  for (df.i in names(df.list)) {
    
    # Loading the DEG tables
    if (length(df.list[[df.i]]) == 0) {
      
      # If there is no object, go to next in loop
      print(paste("There is no", celltype.x, df.i))
      next
      
    } else if (length(df.list[[df.i]]) == 1) {
      
      # Loading the DEG table after parameter
      print(paste("Loading", celltype.x, df.i))
      df = read.xlsx(xlsxFile = df.list[[df.i]])
    }
    
    
    # Adding colnames and coldata
    colnames(df)[1] = "gene"
    df$regulation = "Upregulated"
    
    # Filtering the tables
    df.filt = df[df$p_val_adj < q.value,]
    df.filt = df.filt[df.filt$avg_log2FC > log2fc.cutoff | df.filt$avg_log2FC < -log2fc.cutoff,]
    
    # Editing regulation column after filtering
    df.filt$regulation[df.filt$avg_log2FC < log2fc.cutoff] = "Downregulated"
    
    # Saving the output
    write.xlsx(df.filt, file = paste0(Input.dir, celltype.x, "/", celltype.x, "_", 
                                      df.i, "_", log2fc.cutoff, "_edited_DEG.xlsx"))
    
    # Appending the output to a list of filtered tables
    filtered.tables[[df.i]] = df.filt
  }
  
  if (x.label %in% names(filtered.tables) & y.label %in% names(filtered.tables)) {
    
    # Print what is overlapping
    print(paste("Overlapping", x.label, y.label))
    
    # Merging dataframes to find overlapping genes
    df.overlapped.xy = merge(filtered.tables[[x.label]], filtered.tables[[y.label]], 
                             by = "gene", suffixes = c(paste0("_", x.label), paste0("_", y.label)))
    
    if(nrow(df.overlapped.xy) > 0) {
      
      # Marking restored genes with restored expression
      df.overlapped.xy$restored = "No"
      df.overlapped.xy$restored[df.overlapped.xy[,7] == "Upregulated" & df.overlapped.xy[,13] == "Downregulated"] = "Yes"
      df.overlapped.xy$restored[df.overlapped.xy[,7] == "Downregulated" & df.overlapped.xy[,13] == "Upregulated"] = "Yes"
      
      # Saving the overlapped tables
      write.xlsx(df.overlapped.xy, file = paste0(Input.dir, celltype.x, "/", celltype.x, "_Overlapped_", 
                                                 x.label, "_", y.label, "_", log2fc.cutoff, "_edited_DEG.xlsx"))
    } else if (nrow(df.overlapped.xy) == 0) {
      print("No overlap")
    }
    
  }
  
  if (x.label %in% names(filtered.tables) & z.label %in% names(filtered.tables)) {
    
    # Print what is overlapping
    print(paste("Overlapping", x.label, z.label))
    
    # Merging dataframes to find overlapping genes
    df.overlapped.xz = merge(filtered.tables[[x.label]], filtered.tables[[z.label]], 
                             by = "gene", suffixes = c(paste0("_", x.label), paste0("_", z.label)))
    
    if(nrow(df.overlapped.xz) > 0) {
      
      # Marking restored genes with restored expression
      df.overlapped.xz$restored = "No"
      df.overlapped.xz$restored[df.overlapped.xz[,7] == "Upregulated" & df.overlapped.xz[,13] == "Downregulated"] = "Yes"
      df.overlapped.xz$restored[df.overlapped.xz[,7] == "Downregulated" & df.overlapped.xz[,13] == "Upregulated"] = "Yes"
      
      # Saving the overlapped tables
      write.xlsx(df.overlapped.xz, file = paste0(Input.dir, celltype.x, "/", celltype.x, "_Overlapped_", 
                                                 x.label, "_", z.label, "_", log2fc.cutoff, "_edited_DEG.xlsx"))
    } else if (nrow(df.overlapped.xz) == 0) {
      print("No overlap")
    }
    
  }
  
  if (y.label %in% names(filtered.tables) & z.label %in% names(filtered.tables)) {
    
    # Print what is overlapping
    print(paste("Overlapping", y.label, z.label))
    
    # Merging dataframes to find overlapping genes
    df.overlapped.yz = merge(filtered.tables[[y.label]], filtered.tables[[z.label]], 
                             by = "gene", suffixes = c(paste0("_", y.label), paste0("_", z.label)))
    
    if(nrow(df.overlapped.yz) > 0) {
      
      # Marking restored genes with restored expression
      df.overlapped.yz$restored = "No"
      df.overlapped.yz$restored[df.overlapped.yz[,7] == "Upregulated" & df.overlapped.yz[,13] == "Downregulated"] = "Yes"
      df.overlapped.yz$restored[df.overlapped.yz[,7] == "Downregulated" & df.overlapped.yz[,13] == "Upregulated"] = "Yes"
      
      # Saving the overlapped tables
      write.xlsx(df.overlapped.yz, file = paste0(Input.dir, celltype.x, "/", celltype.x, "_Overlapped_", 
                                                 y.label, "_", z.label, "_", log2fc.cutoff, "_edited_DEG.xlsx"))
    } else if (nrow(df.overlapped.yz) == 0) {
      print("No overlap")
    }
    
  }
  
  # Generate a Venndiagram of overlapping genes
  xyz.genelists = list(x.genelist = filtered.tables[[x.label]][[1]],
                       y.genelist = filtered.tables[[y.label]][[1]],
                       z.genelist = filtered.tables[[z.label]][[1]])
  
  names(xyz.genelists) = c(x.label, y.label, z.label)
  
  vennplot = ggvenn(xyz.genelists,
                    #fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"), 
                    fill_color = Venn.cols,
                    stroke_size = 0.5, set_name_size = 5, text_size = 3)
  
  ggsave2(plot = vennplot, filename = paste0(Input.dir, celltype.x, "/", celltype.x,"_Venndiagram_DEGs.pdf"), dpi = 700)
  
}


GO_enrichment_analysis <- function(celltype.x, skip.overlap = FALSE) {
  
  # p and q values for GO enrichment analysis
  p = 0.05
  q = 0.05
  p.KEGG = 0.1 
  
  # Setting input and output
  input.dir.list = list.files(path = paste0(Input.dir, celltype.x), pattern = "*_edited_DEG.xlsx", full.names = TRUE)
  DEG.comp = list.files(path = paste0(Input.dir, celltype.x), pattern = "*_edited_DEG.xlsx", full.names = FALSE)
  DEG.comp = str_remove(DEG.comp, pattern = ".xlsx")

  # Read the DEG table
  DEG.df.list = lapply(input.dir.list, read.xlsx)
  DEG.comp = sub("([A-Za-z]+_[A-Za-z]+)_0.5.*$", "\\1", DEG.comp)
  names(DEG.df.list) = DEG.comp
  
  print(paste0("GO enrichment analysis on ", celltype.x))
  
  # Checking and/or generating output dir
  if (dir.exists(path = paste0(Input.dir, celltype.x, "/GSEA_enrichment_analysis/")) == FALSE) {
    print(paste0("Generating output directory ", Input.dir, celltype.x, "/GSEA_enrichment_analysis/"))
    dir.create(path = paste0(Input.dir, celltype.x, "/GSEA_enrichment_analysis/"), recursive = TRUE)
    Output.dir.GSEA = paste0(Input.dir, celltype.x, "/GSEA_enrichment_analysis/")
  } else if (dir.exists(path = paste0(Input.dir, celltype.x, "/GSEA_enrichment_analysis/")) == TRUE) {
    print("Directory exists")
    Output.dir.GSEA = paste0(Input.dir, celltype.x, "/GSEA_enrichment_analysis/")
  } else {
    print("Error with output directory")
  }
  
  # Running GO on each dataframe
  for (x in DEG.comp) {
    
    # Flag if the input is overlapped or not
    Overlap.flag = grepl("_Overlapped_", x, fixed = TRUE)
    
    if (skip.overlap == TRUE & Overlap.flag == TRUE) {
      print("Skipping overlap df")
      next
    }

    # Generate lists to append enrichGO and compareCluster output to
    # for later plotting
    enrichGO.list = list()
    enrichGO.names = c()
    compareCluster.list = list()
    compareCluster.names = c()
    
    print(paste0("GO enrichment analysis on ", x))
    
    if (is.null(DEG.df.list[[x]])) {
      print("Skipping the df, it is NULL")
      next
      
    }
    
    df = DEG.df.list[[x]]
    
    if (nrow(df) < 5) {
      print("Skipping the df, less than 5 rows")
      next
    }
    
    colnames(df)[1] = "SYMBOL"
    colnames(df)[3] = "avg_log2FC"

    if (Overlap.flag == FALSE) {
      colnames(df)[7] = "group"
    } else if (Overlap.flag == TRUE) {
      colnames(df)[14] = "group"
    }
    
    # Split into down and upregulated genes
    df.down = df[df$group == "Downregulated",]
    df.up = df[df$group == "Upregulated",]
    
    # Setting the output dir
    Output.dir = paste0(Output.dir.GSEA, x)
    
    # Run enrichGO from clusterprofiler on on biological processes, molecular function
    # ceullular compartments and all
    
    #### Biological processes #######
    enrichGO.BP.all = enrichGO(gene = df$SYMBOL, 
                               OrgDb = org.Hs.eg.db,
                               keyType = 'SYMBOL',
                               ont = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = p,
                               qvalueCutoff  = q)
    
    if (is.null(x = enrichGO.BP.all) == FALSE) {
      
      #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
      enrichGO.BP.all@result = na.omit(enrichGO.BP.all@result)
      
      if (nrow(enrichGO.BP.all@result) > 1) {
        
        df.output = enrichGO.BP.all@result
        df.output$enriched = "No"
        df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichGO.BP.all_table.xlsx"))
        
        enrichGO.list = append(enrichGO.list, enrichGO.BP.all, after = length(enrichGO.list))
        enrichGO.names[length(enrichGO.list)] = "BP"
         
      }
      
    }
    
    if (Overlap.flag == FALSE & skip.overlap == FALSE) {
      enrichGO.BP.all.compare = compareCluster(SYMBOL~group, data=df, fun="enrichGO",
                                               OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
                                               ont = "BP", pAdjustMethod = "BH",
                                               pvalueCutoff  = p, qvalueCutoff  = q)
      
      if (is.null(x = enrichGO.BP.all.compare) == FALSE) {
        
        #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
        enrichGO.BP.all.compare@compareClusterResult = na.omit(enrichGO.BP.all.compare@compareClusterResult)
        
        if (nrow(enrichGO.BP.all.compare@compareClusterResult) > 1) {
          df.output = enrichGO.BP.all.compare@compareClusterResult
          df.output$enriched = "No"
          df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
          openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_compareCluster.BP.all_table.xlsx"))
          
          compareCluster.list = append(compareCluster.list, enrichGO.BP.all.compare, after = length(compareCluster.list))
          compareCluster.names[length(compareCluster.list)] = "BP"
        }
        
      }
    }
    
    #### Molecular functions #######
    enrichGO.MF.all = enrichGO(gene = df$SYMBOL, 
                               OrgDb = org.Hs.eg.db,
                               keyType = 'SYMBOL',
                               ont = "MF",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = p,
                               qvalueCutoff  = q)
    
    if (is.null(x =enrichGO.MF.all) == FALSE & nrow(enrichGO.MF.all@result) > 1) {
      
      #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
      enrichGO.MF.all@result = na.omit(enrichGO.MF.all@result)
      
      if (nrow(enrichGO.MF.all@result) > 1) {
        
        df.output = enrichGO.MF.all@result
        df.output$enriched = "No"
        df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichGO.MF.all_table.xlsx"))
        
        enrichGO.list = append(enrichGO.list, enrichGO.MF.all, after = length(enrichGO.list))
        enrichGO.names[length(enrichGO.list)] = "MF"
        
      }
      
    }
    
    if (Overlap.flag == FALSE & skip.overlap == FALSE) {
      enrichGO.MF.all.compare = compareCluster(SYMBOL~group, data=df, fun="enrichGO",
                                               OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
                                               ont = "MF", pAdjustMethod = "BH",
                                               pvalueCutoff  = p, qvalueCutoff  = q)
      
      if (is.null(x = enrichGO.MF.all.compare) == FALSE) {
        
        #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
        enrichGO.MF.all.compare@compareClusterResult = na.omit(enrichGO.MF.all.compare@compareClusterResult)
        
        if (nrow(enrichGO.MF.all.compare@compareClusterResult) > 1) {
          df.output = enrichGO.MF.all.compare@compareClusterResult
          df.output$enriched = "No"
          df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
          openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_compareCluster.MF.all_table.xlsx"))
          
          compareCluster.list = append(compareCluster.list, enrichGO.MF.all.compare, after = length(compareCluster.list))
          compareCluster.names[length(compareCluster.list)] = "MF"
        }

      }
      
    }
    
    #### Cellular compartments #######
    enrichGO.CC.all = enrichGO(gene = df$SYMBOL, 
                               OrgDb = org.Hs.eg.db,
                               keyType = 'SYMBOL',
                               ont = "CC",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = p,
                               qvalueCutoff  = q)
    
    if (is.null(x =enrichGO.CC.all) == FALSE) {
      
      #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
      enrichGO.CC.all@result = na.omit(enrichGO.CC.all@result)
      
      if (nrow(enrichGO.CC.all@result) > 1) {
        df.output = enrichGO.CC.all@result
        df.output$enriched = "No"
        df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichGO.CC.all_table.xlsx"))
        
        enrichGO.list = append(enrichGO.list, enrichGO.CC.all, after = length(enrichGO.list))
        enrichGO.names[length(enrichGO.list)] = "CC"
      }
      
    }
    
    if (Overlap.flag == FALSE & skip.overlap == FALSE) {
      enrichGO.CC.all.compare = compareCluster(SYMBOL~group, data=df, fun="enrichGO",
                                               OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
                                               ont = "CC", pAdjustMethod = "BH",
                                               pvalueCutoff  = p, qvalueCutoff  = q)
      
      if (is.null(x = enrichGO.CC.all.compare) == FALSE) {
        
        #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
        enrichGO.CC.all.compare@compareClusterResult = na.omit(enrichGO.CC.all.compare@compareClusterResult)
        
        if (nrow(enrichGO.CC.all.compare@compareClusterResult) > 1) {
          df.output = enrichGO.CC.all.compare@compareClusterResult
          df.output$enriched = "No"
          df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
          openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_compareCluster.CC.all_table.xlsx"))
          
          compareCluster.list = append(compareCluster.list, enrichGO.CC.all.compare, after = length(compareCluster.list))
          compareCluster.names[length(compareCluster.list)] = "CC"
          
        }
        
      }
    }
    
    #### All subontologies #######
    enrichGO.ALL.all = enrichGO(gene = df$SYMBOL, 
                               OrgDb = org.Hs.eg.db,
                               keyType = 'SYMBOL',
                               ont = "ALL",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = p,
                               qvalueCutoff  = q)
    
    if (is.null(x = enrichGO.ALL.all) == FALSE) {
      
      #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
      enrichGO.ALL.all@result = na.omit(enrichGO.ALL.all@result)
      
      if (nrow(enrichGO.ALL.all@result) > 1) {
        df.output = enrichGO.ALL.all@result
        df.output$enriched = "No"
        df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichGO.ALL.all_table.xlsx"))
        
        enrichGO.list = append(enrichGO.list, enrichGO.ALL.all, after = length(enrichGO.list))
        enrichGO.names[length(enrichGO.list)] = "ALL" 
      }
      
    }
    
    if (Overlap.flag == FALSE & skip.overlap == FALSE) {
      enrichGO.ALL.all.compare = compareCluster(SYMBOL~group, data=df, fun="enrichGO",
                                                OrgDb = org.Hs.eg.db, keyType = 'SYMBOL',
                                                ont = "ALL", pAdjustMethod = "BH",
                                                pvalueCutoff  = p, qvalueCutoff  = q)
      
      if (is.null(x = enrichGO.ALL.all.compare) == FALSE) {
        
        #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
        enrichGO.ALL.all.compare@compareClusterResult = na.omit(enrichGO.ALL.all.compare@compareClusterResult)
        
        if (nrow(enrichGO.ALL.all.compare@compareClusterResult) > 1) {
    
          df.output = enrichGO.ALL.all.compare@compareClusterResult
          df.output$enriched = "No"
          df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
          openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_compareCluster.ALL.all_table.xlsx"))
          
          compareCluster.list = append(compareCluster.list, enrichGO.ALL.all.compare, after = length(compareCluster.list))
          compareCluster.names[length(compareCluster.list)] = "ALL"
        }
        
      }
      
    }
    
    #### GO Hallmarks ######
    enrichGO.Hallmark.all = enricher(gene = df$SYMBOL, 
                                     TERM2GENE= m_Hallmark,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = p,
                                     qvalueCutoff  = q)
    
    if (is.null(x = enrichGO.Hallmark.all) == FALSE) {
      
      #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
      enrichGO.Hallmark.all@result = na.omit(enrichGO.Hallmark.all@result)
      
      if (nrow(enrichGO.Hallmark.all@result) > 1) {
        df.output = enrichGO.Hallmark.all@result
        df.output$enriched = "No"
        df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichGO.Hallmark.all_table.xlsx"))
        
        enrichGO.list = append(enrichGO.list, enrichGO.Hallmark.all, after = length(enrichGO.list))
        enrichGO.names[length(enrichGO.list)] = "Hallmark" 
      }
      
    }
    
    if (Overlap.flag == FALSE & skip.overlap == FALSE) {
      enrichGO.Hallmark.all.compare = compareCluster(SYMBOL~group, data=df, fun="enricher",
                                                     TERM2GENE= m_Hallmark, pAdjustMethod = "BH",
                                                     pvalueCutoff  = p, qvalueCutoff  = q)
      
      if (is.null(x = enrichGO.Hallmark.all.compare) == FALSE) {
        
        #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
        enrichGO.Hallmark.all.compare@compareClusterResult = na.omit(enrichGO.Hallmark.all.compare@compareClusterResult)
        
        if (nrow(enrichGO.Hallmark.all.compare@compareClusterResult) > 1) {
          df.output = enrichGO.Hallmark.all.compare@compareClusterResult
          df.output$enriched = "No"
          df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
          openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_compareCluster.Hallmark.all_table.xlsx"))
          
          compareCluster.list = append(compareCluster.list, enrichGO.Hallmark.all.compare, after = length(compareCluster.list))
          compareCluster.names[length(compareCluster.list)] = "Hallmark"
        }
        
      }
    }
    
    #### GO Immune ######
    enrichGO.Immune.all = enricher(gene = df$SYMBOL, 
                                     TERM2GENE= m_Immune,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = p,
                                     qvalueCutoff  = q)
    
    if (is.null(x = enrichGO.Immune.all) == FALSE) {
      
      #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
      enrichGO.Immune.all@result = na.omit(enrichGO.Immune.all@result)
      
      if (nrow(enrichGO.Immune.all@result) > 1) {
        df.output = enrichGO.Immune.all@result
        df.output$enriched = "No"
        df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichGO.Immune.all_table.xlsx"))
        
        enrichGO.list = append(enrichGO.list, enrichGO.Immune.all, after = length(enrichGO.list))
        enrichGO.names[length(enrichGO.list)] = "Immune" 
      }
      
    }
    
    if (Overlap.flag == FALSE & skip.overlap == FALSE) {
      enrichGO.Immune.all.compare = compareCluster(SYMBOL~group, data=df, fun="enricher",
                                                     TERM2GENE= m_Immune, pAdjustMethod = "BH",
                                                     pvalueCutoff  = p, qvalueCutoff  = q)
      
      if (is.null(x = enrichGO.Immune.all.compare) == FALSE) {
        
        #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
        enrichGO.Immune.all.compare@compareClusterResult = na.omit(enrichGO.Immune.all.compare@compareClusterResult)
        
        if (nrow(enrichGO.Immune.all.compare@compareClusterResult) > 1) {
          df.output = enrichGO.Immune.all.compare@compareClusterResult
          df.output$enriched = "No"
          df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
          openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_compareCluster.Immune.all_table.xlsx"))
          
          compareCluster.list = append(compareCluster.list, enrichGO.Immune.all.compare, after = length(compareCluster.list))
          compareCluster.names[length(compareCluster.list)] = "Immune"
        }

      }
    }
    
    #### GO HPO ######
    enrichGO.HPO.all = enricher(gene = df$SYMBOL, 
                                   TERM2GENE= m_HPO,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = p,
                                   qvalueCutoff  = q)
    
    if (is.null(x = enrichGO.HPO.all) == FALSE) {
      
      #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
      enrichGO.HPO.all@result = na.omit(enrichGO.HPO.all@result)
      
      if (nrow(enrichGO.HPO.all@result) > 1) {
        df.output = enrichGO.HPO.all@result
        df.output$enriched = "No"
        df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichGO.HPO.all_table.xlsx"))
        
        enrichGO.list = append(enrichGO.list, enrichGO.HPO.all, after = length(enrichGO.list))
        enrichGO.names[length(enrichGO.list)] = "HPO" 
      }
      
    }
    
    if (Overlap.flag == FALSE & skip.overlap == FALSE) {
      enrichGO.HPO.all.compare = compareCluster(SYMBOL~group, data=df, fun="enricher",
                                                   TERM2GENE= m_HPO, pAdjustMethod = "BH",
                                                   pvalueCutoff  = p, qvalueCutoff  = q)
      
      if (is.null(x = enrichGO.HPO.all.compare) == FALSE) {
        
        #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
        enrichGO.HPO.all.compare@compareClusterResult = na.omit(enrichGO.HPO.all.compare@compareClusterResult)
        
        if (nrow(enrichGO.HPO.all.compare@compareClusterResult) > 1) {
          df.output = enrichGO.HPO.all.compare@compareClusterResult
          df.output$enriched = "No"
          df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
          openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_compareCluster.HPO.all_table.xlsx"))
          
          compareCluster.list = append(compareCluster.list, enrichGO.HPO.all.compare, after = length(compareCluster.list))
          compareCluster.names[length(compareCluster.list)] = "HPO"
        }
        
      }
    }
    
    #### GO msigdb_BP ######
    enrichGO.msigdb_BP.all = enricher(gene = df$SYMBOL, 
                                TERM2GENE= m_msigdb_BP,
                                pAdjustMethod = "BH",
                                pvalueCutoff  = p,
                                qvalueCutoff  = q)
    
    if (is.null(x = enrichGO.msigdb_BP.all) == FALSE) {
      
      #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
      enrichGO.msigdb_BP.all@result = na.omit(enrichGO.msigdb_BP.all@result)
      
      if (nrow(enrichGO.msigdb_BP.all@result) > 1) {
        df.output = enrichGO.msigdb_BP.all@result
        df.output$enriched = "No"
        df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichGO.msigdb_BP.all_table.xlsx"))
        
        enrichGO.list = append(enrichGO.list, enrichGO.msigdb_BP.all, after = length(enrichGO.list))
        enrichGO.names[length(enrichGO.list)] = "msigdb_BP" 
      }
      
    }
    
    if (Overlap.flag == FALSE & skip.overlap == FALSE) {
      enrichGO.msigdb_BP.all.compare = compareCluster(SYMBOL~group, data=df, fun="enricher",
                                                TERM2GENE= m_msigdb_BP, pAdjustMethod = "BH",
                                                pvalueCutoff  = p, qvalueCutoff  = q)
      
      if (is.null(x = enrichGO.msigdb_BP.all.compare) == FALSE) {
        
        #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
        enrichGO.msigdb_BP.all.compare@compareClusterResult = na.omit(enrichGO.msigdb_BP.all.compare@compareClusterResult)
        
        if (nrow(enrichGO.msigdb_BP.all.compare@compareClusterResult) > 1) {
          df.output = enrichGO.msigdb_BP.all.compare@compareClusterResult
          df.output$enriched = "No"
          df.output$enriched[df.output$p.adjust <= p & df.output$qvalue <= q] = "Yes"
          openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_compareCluster.msigdb_BP.all_table.xlsx"))
          
          compareCluster.list = append(compareCluster.list, enrichGO.msigdb_BP.all.compare, after = length(compareCluster.list))
          compareCluster.names[length(compareCluster.list)] = "msigdb_BP"
        }
        
      }
    }
    
    ##### Run KEGG pathway over-representation analysis ####
    df$ENTREZID = mapIds(org.Hs.eg.db, keys = df$SYMBOL, keytype = "SYMBOL", column = "ENTREZID")
    
    enrichKEGG.all <- enrichKEGG(gene = df$ENTREZID,
                     organism = 'hsa',
                     pAdjustMethod = "BH",
                     pvalueCutoff = p.KEGG,
                     qvalueCutoff = q,
                     keyType = "kegg")
    
    #if (is.null(x = enrichKEGG.all) == FALSE & nrow(enrichKEGG.all@result) > 1) {
    if (is.null(x = enrichKEGG.all) == FALSE) {
      
      #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
      enrichKEGG.all <- setReadable(enrichKEGG.all, 'org.Hs.eg.db', 'ENTREZID')
      enrichKEGG.all@result = na.omit(enrichKEGG.all@result)
      
      if (nrow(enrichKEGG.all@result) > 1) {
        df.output = enrichKEGG.all@result
        df.output$enriched = "No"
        df.output$enriched[df.output$p.adjust <= p.KEGG & df.output$qvalue <= q] = "Yes"
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichKEGG.all_table.xlsx")) 
      }
      
    }
    
    # Split the object in up and downregulated dataframes to run seperate analysis
    if (Overlap.flag == FALSE & skip.overlap == FALSE) {
      
      df.down = df[df$group == "Downregulated",]
      df.up = df[df$group == "Upregulated",]
      
      if (nrow(df.down) > 1) {
        enrichKEGG.down <- enrichKEGG(gene = df.down$ENTREZID,
                                     organism = 'hsa',
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = p.KEGG,
                                     qvalueCutoff = q,
                                     keyType = "kegg")
        
        if (is.null(x = enrichKEGG.down) == FALSE) {
          
          #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
          enrichKEGG.down <- setReadable(enrichKEGG.down, 'org.Hs.eg.db', 'ENTREZID')
          enrichKEGG.down@result = na.omit(enrichKEGG.down@result)
          
          if (nrow(enrichKEGG.down@result) > 1) {
            df.output.down = enrichKEGG.down@result
            df.output.down$enriched = "No"
            df.output.down$enriched[df.output.down$p.adjust <= p.KEGG & df.output.down$qvalue <= q] = "Yes"
            df.output.down$group = "Downregulated"
            openxlsx::write.xlsx(df.output.down, file = paste0(Output.dir, "_enrichKEGG.down_table.xlsx"))
          }
          
          
        } else if (is.null(x = enrichKEGG.down) == TRUE | nrow(enrichKEGG.down@result) == 0) {
          df.output.down = NULL
        }
      }
        
      if (nrow(df.up) > 1) {
        
        enrichKEGG.up <- enrichKEGG(gene = df.up$ENTREZID,
                                      organism = 'hsa',
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = p.KEGG,
                                      qvalueCutoff = q,
                                      keyType = "kegg")
        
        if (is.null(x = enrichKEGG.up) == FALSE) {
          
          #Remove rows in output with NA values, often caused by no q-value i.e. weak predictions
          enrichKEGG.up <- setReadable(enrichKEGG.up, 'org.Hs.eg.db', 'ENTREZID')
          enrichKEGG.up@result = na.omit(enrichKEGG.up@result)
          
          if (nrow(enrichKEGG.up@result) > 1) {
            df.output.up = enrichKEGG.up@result
            df.output.up$enriched = "No"
            df.output.up$enriched[df.output.up$p.adjust <= p.KEGG & df.output.up$qvalue <= q] = "Yes"
            df.output.up$group = "Upregulated" 
          }
          
        } else if (is.null(x = enrichKEGG.up) == TRUE | nrow(enrichKEGG.up@result) == 0) {
          df.output.up = NULL
        }
        
      }
      
      if (!is.null(df.output.down) & !is.null(df.output.up)) {
        df.output = rbind(df.output.down, df.output.up)
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichKEGG.all_regulated_table.xlsx"))
      } else if (!is.null(df.output.down) & is.null(df.output.up)) {
        df.output = df.output.down
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichKEGG.all_regulated_table.xlsx"))
      } else if (is.null(df.output.down) & !is.null(df.output.up)) {
        df.output = df.output.up
        openxlsx::write.xlsx(df.output, file = paste0(Output.dir, "_enrichKEGG.all_regulated_table.xlsx"))
      }
      
    }

    # Add names to generated lists
    print(enrichGO.names)
    print(compareCluster.names)
    
    names(enrichGO.list) = enrichGO.names
    names(compareCluster.list) = compareCluster.names
    
    # Loop over the enrichGO list. Skip if there are no significant enrichments
    for (i in enrichGO.names) {
      print(i)
      enrich.df = enrichGO.list[[i]]
     
     if (sum(enrich.df@result$qvalue < q) >= 1 & sum(enrich.df@result$p.adjust < p) >= 1) {
       
       print(paste("Plotting", i))
       
       plot.res = barplot(enrich.df, title = paste0("GO ", i, " in ", celltype.x), showCategory = 10) + xlim(0,10) + 
         theme_cowplot()
       ggsave2(plot = plot.res, filename = paste0(Output.dir, "_enrichGO.", i,".all_barplot.pdf"), dpi = 700)
       
       log.enrinch.df = enrich.df@result
       log.enrinch.df$log10_padjust = -log10(log.enrinch.df$p.adjust)
       log.enrinch.df = log.enrinch.df[order(log.enrinch.df$log10_padjust, decreasing = TRUE),]
       
       plot.res = ggplot(log.enrinch.df[1:10,], aes(x = reorder(Description, log10_padjust), y = log10_padjust)) + 
         geom_bar(stat = "identity", fill = "red") + ylim(0, 4) + coord_flip() + theme_cowplot()
       ggsave2(plot = plot.res, filename = paste0(Output.dir, "_enrichGO.", i,".all_log10_padjust_barplot.pdf"), 
               dpi = 700, width = 10, height = 3)
       
       plot.res = dotplot(enrich.df, title = paste0("GO ", i, " in ", celltype.x), showCategory = 10) + theme_cowplot()
       ggsave2(plot = plot.res, filename = paste0(Output.dir, "_enrichGO.", i,".all_dotplot.pdf"), dpi = 700)
       
       plot.res = heatplot(enrich.df, foldChange=df$avg_log2FC, showCategory=10) + theme_cowplot() +
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "italic")) + 
         labs(title = paste0("GO ", i, " in ", celltype.x))
       ggsave2(plot = plot.res, filename = paste0(Output.dir, "enrichGO.", i,".all_heatplot.pdf"),
               width = 25, height = 8, dpi = 700)
       
     } else if (sum(enrich.df@result$qvalue < q) < 1 | sum(enrich.df@result$p.adjust < p) < 1) {
       print("No enriched results")
     }
     
    }
    
    # Loop over the compareCLuster list. Skip if there are no significant enrichments
    if (Overlap.flag == FALSE & skip.overlap == FALSE) {
      for (i in compareCluster.names) {
        print(i)
        enrich.df = compareCluster.list[[i]]
        
       if (sum(enrich.df@compareClusterResult$qvalue < q) >= 1 & 
            sum(enrich.df@compareClusterResult$p.adjust < p) >= 1) {
          
          print(paste("Plotting", i))
          
          plot.res = dotplot(enrich.df, showCategory = 10) + theme_cowplot() +
            labs(title = paste0("GO ", i, " in ", celltype.x))
          ggsave2(plot = plot.res, filename = paste0(Output.dir, "_compareCluster.", i,".dotplot.pdf"),
                  height = 7 , width = 8, dpi = 700)
          
          # For cnetplot and emapplot, there have to be more than one rows in enrich.df results
          if (sum(enrich.df@compareClusterResult$qvalue < q) > 1 & 
              sum(enrich.df@compareClusterResult$p.adjust < p) > 1) {
            
            plot.res = cnetplot(enrich.df, showCategory = 10, node_label = "all") +
              scale_fill_manual(values=Two.cols)
            ggsave2(plot = plot.res, filename = paste0(Output.dir, "_compareCluster.", i,".cnetplot.pdf"),dpi = 700)
            
            enrich.df.pair <- pairwise_termsim(enrich.df)
            plot.res = emapplot(enrich.df.pair, showCategory = 20) + scale_fill_manual(values=Two.cols)
            ggsave2(plot = plot.res, filename = paste0(Output.dir, "_compareCluster.", i,".emapplot.pdf"),dpi = 700,
                    width = 12, height = 12)
            
          }
          
        } else if (sum(enrich.df@compareClusterResult$qvalue < q) < 1 | 
                   sum(enrich.df@compareClusterResult$p.adjust < p) < 1) {
          print("No enriched results")
        }
    }
      
      print(paste0("Done with ", i, " in ", celltype.x))
      
    }
    
  }
  
}

DEG.GO.organiser <- function(celltype.x) {
  
  # Setting input and output
  input.dir.list = list.files(path = paste0(Input.dir, celltype.x), pattern = "*.all_table.xlsx", full.names = TRUE)
  DEG.comp = list.files(path = paste0(Input.dir, celltype.x), pattern = "*.all_table.xlsx", full.names = FALSE)
  
  # Read the DEG table
  DEG.df.list = lapply(input.dir.list, read.xlsx)
  DEG.comp = sub("([A-Za-z]+_[A-Za-z]+)_0.5.*$", "\\1", DEG.comp)
  names(DEG.df.list) = DEG.comp
  
  print(paste("Organising the DEG and GO lists in", celltype.x))
  
  # Running GO on each dataframe
  for (x in DEG.comp) {
    
    print(paste("Organising", x))
    
  }
  
}

DEG.overlap.plotting <- function(seurat.x, celltype.x, selected.genes = "ESR1", group.name = "all",
                                 treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                 scale.min.x = 20, scale.max.x = 100, dot.scale.x = 8) {
  
  # Subsetting target celltype.x
  Idents(object = seurat.x) <- select.idents
  subset.cell = subset(seurat.x, idents = celltype.x)
  
  # Normalizing the data
  DefaultAssay(subset.cell) = "RNA"
  Idents(subset.cell) = "orig.ident"
  subset.cell = NormalizeData(subset.cell)
  subset.cell = FindVariableFeatures(subset.cell)
  all.genes = rownames(subset.cell)
  subset.cell = ScaleData(subset.cell, features = all.genes)  #endo.integrated = FindVariableFeatures(endo.integrated)
  Idents(subset.cell) = "Group_Stage"
  DefaultAssay(subset.cell) = "RNA"
  
  # Saving the subset object if true
  if (save.subset.seurat == TRUE) {
    SaveH5Seurat(subset.cell, paste0(Output.dir,Project_name, "_", celltype.x, "_seurat.h5seurat"),
                 overwrite = TRUE)
  }
  
  # Re-order the object for plotting. Make it Control vs. PCOS in Alt 2
  
  subset.cell$Group_Stage <- factor(subset.cell$Group_Stage, levels = treatment.levels)
  Idents(subset.cell) <- "Group_Stage"
  
  subset.cell.alt2 = subset.cell
  subset.cell.alt2$Group_Stage <- factor(subset.cell.alt2$Group_Stage, levels = rev(treatment.levels))
  Idents(subset.cell.alt2) <- "Group_Stage"
  
  # Generating dotplots
  plot.res = DotPlot(subset.cell, features = selected.genes, scale = FALSE, col.min = 0, col.max = 5, scale.min = scale.min.x, scale.max = scale.max.x,
                     cols = Two.cols, dot.scale = dot.scale.x, assay = "RNA") + theme(axis.text.y = element_text(face = "italic")) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) + RotatedAxis() + coord_flip()
  ggsave2(plot = plot.res, filename = paste0(Output.dir,Project_name, "_", celltype.x, "_", group.name, "_selected_genes_Alt1_Dotplot.pdf"), dpi = 700)
  
  plot.res = DotPlot(subset.cell.alt2, features = selected.genes, scale = FALSE, col.min = 0, col.max = 5, scale.min = scale.min.x, scale.max = scale.max.x,
                     cols = Two.cols, dot.scale = dot.scale.x, assay = "RNA") + theme(axis.text.x = element_text(face = "italic")) +
    theme(axis.text.x=element_text(angle = 45, hjust = 1)) 
  ggsave2(plot = plot.res, filename = paste0(Output.dir,Project_name, "_", celltype.x, "_", group.name, "_selected_genes_Alt2_Dotplot.pdf"), dpi = 700)
  
  
  # Generating violin plots
  plot.res = VlnPlot(subset.cell, features = selected.genes, pt.size = 0, cols = Group.cols) & 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          plot.title = element_text(face = "italic"))
  ggsave2(plot = plot.res, paste0(Output.dir,Project_name, "_", celltype.x, "_", group.name, "_selected_genes_Vlnplot.pdf"), dpi = 700)
  #ggsave2(plot = plot.res, paste0(Output.dir.,Project_name, "_", celltype.x, "_selected_genes_Vlnplot.pdf"), dpi = 700, 
  #        height = 10, width = 10)
  
}

# Subsetting and plotting the seurat object is required to be done on the
# UPPMAX server
if (server.run == TRUE) {
  
  # Loading reclustered Seurat
  print("Seurat object loading")
  endo.integrated = LoadH5Seurat(file = paste0(seurat.dir,seurat.project, "_reclustered_labelled.h5seurat"))
  Idents(object = endo.integrated) <- "Group_Stage"
  endo.integrated.Ctrl_PCOS = subset(endo.integrated, idents = c("Control", "PCOS_W0"))

  if (select.idents == "Epithelium_labelled") {
    # Plotting control vs. PCOS
    lumenal.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                celltype.x = "Lumenal", group.name = "Ctrl_PCOS",
                                                treatment.levels = c("Control", "PCOS_W0"),
                                                selected.genes = c("SEMA3E", "ROBO2", "FGF13", "KIAA1217", "CTNNA2", #Top 5 downregulated
                                                                   "INPP4B", "RIMKLB", "CPM", "SLPI", "PAEP")) # Top 5 upregulated
    
    ciliated.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                 celltype.x = "Ciliated", group.name = "Ctrl_PCOS",
                                                 treatment.levels = c("Control", "PCOS_W0"),
                                                 selected.genes = c("CCSER1", "CERNA2", "MS4A8", "RNLS", "RP1", #Top 5 downregulated
                                                                    "PLCE1", "LCN2", "CCDC171", "CCNO", "CDC20B")) # Top 5 upregulated
    
    AR.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                           celltype.x = "AR+", group.name = "Ctrl_PCOS",
                                           treatment.levels = c("Control", "PCOS_W0"),
                                           selected.genes = c("ROBO2", "SFRP1", "SEMA3E", "LINC00621", "LINC02456", #Top 5 downregulated
                                                              "STC1", "NEAT1", "PAEP", "MMP10", "MMP3")) # Top 5 upregulated
    
    SOX9prolif.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                   celltype.x = "SOX9+ proliferative", group.name = "Ctrl_PCOS",
                                                   treatment.levels = c("Control", "PCOS_W0"),
                                                   selected.genes = c("SEMA3E", "ROBO2", "HMCN1", "SOX5", "LGR5", #Top 5 downregulated
                                                                      "SLC1A1", "A2ML1-AS1", "CATSPERB", "RIMKLB", "PAEP")) # Top 5 upregulated
    
    SOX9posLGR5neg.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                       celltype.x = "SOX9+ LGR5-", group.name = "Ctrl_PCOS",
                                                       treatment.levels = c("Control", "PCOS_W0"),
                                                       selected.genes = c("SEMA3E", "GPC6", "LGR5", "CTNNA2", "NRCAM", #Top 5 downregulated
                                                                          "NEAT1", "ARHGAP26", "RIMKLB", "SLPI", "PAEP")) # Top 5 upregulated
    
    SOX9posLGR5pos.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                       celltype.x = "SOX9+ LGR5+", group.name = "Ctrl_PCOS",
                                                       treatment.levels = c("Control", "PCOS_W0"),
                                                       selected.genes = c("HMCN1", "SEMA3E", "ROBO2", "KIAA1217", "DIAPH3", #Top 5 downregulated
                                                                          "SLPI", "ITGA2", "LAMC2", "ADAMTS9", "PAEP")) # Top 5 upregulated
    
    # Plotting all groups
    lumenal.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                celltype.x = "Lumenal", group.name = "Top10",
                                                treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                                selected.genes = c("SEMA3E", "ROBO2", "FGF13", "KIAA1217", "CTNNA2", #Top 5 downregulated
                                                                   "INPP4B", "RIMKLB", "CPM", "SLPI", "PAEP")) # Top 5 upregulated
  
    ciliated.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                 celltype.x = "Ciliated", group.name = "Top10",
                                                 treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                                 selected.genes = c("CCSER1", "CERNA2", "MS4A8", "RNLS", "RP1", #Top 5 downregulated
                                                                    "PLCE1", "LCN2", "CCDC171", "CCNO", "CDC20B")) # Top 5 upregulated
    
    AR.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                           celltype.x = "AR+", group.name = "Top10",
                                           treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                           selected.genes = c("ROBO2", "SFRP1", "SEMA3E", "LINC00621", "LINC02456", #Top 5 downregulated
                                                              "STC1", "NEAT1", "PAEP", "MMP10", "MMP3")) # Top 5 upregulated
    
    SOX9prolif.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                   celltype.x = "SOX9+ proliferative", group.name = "Top10",
                                                   treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                                   selected.genes = c("SEMA3E", "ROBO2", "HMCN1", "SOX5", "LGR5", #Top 5 downregulated
                                                                      "SLC1A1", "A2ML1-AS1", "CATSPERB", "RIMKLB", "PAEP")) # Top 5 upregulated
    
    SOX9posLGR5neg.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                       celltype.x = "SOX9+ LGR5-", group.name = "Top10",
                                                       treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                                       selected.genes = c("SEMA3E", "GPC6", "LGR5", "CTNNA2", "NRCAM", #Top 5 downregulated
                                                                          "NEAT1", "ARHGAP26", "RIMKLB", "SLPI", "PAEP")) # Top 5 upregulated
    
    SOX9posLGR5pos.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                       celltype.x = "SOX9+ LGR5+", group.name = "Top10",
                                                       treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                                       selected.genes = c("HMCN1", "SEMA3E", "ROBO2", "KIAA1217", "DIAPH3", #Top 5 downregulated
                                                                          "SLPI", "ITGA2", "LAMC2", "ADAMTS9", "PAEP")) # Top 5 upregulated
    
    # Plotting restored genes PCOS Metformin
    lumenal.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                celltype.x = "Lumenal", group.name = "Restored",
                                                treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                                selected.genes = c("GREM2", "GRAMD1C", "NAV3", "IDO1", "ADGRL2", #Top 5 downregulated
                                                                   "RGS6", "INPP4B", "RIMKLB", "SLPI", "PAEP")) # Top 5 upregulated
    
    ciliated.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                 celltype.x = "Ciliated", group.name = "Restored",
                                                 treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                                 selected.genes = c("LRMP", "NEAT1", "SOD2")) # Top 5 upregulated
    
    AR.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                           celltype.x = "AR+", group.name = "Restored",
                                           treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                           selected.genes = c("ROBO2", "CTNNA2", "FILIP1L", "MIR924HG", "SFRP4", #Top 5 downregulated
                                                              "STC1", "NEAT1", "PAEP", "MMP10", "MMP3")) # Top 5 upregulated
    
    SOX9prolif.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                   celltype.x = "SOX9+ proliferative", group.name = "Restored",
                                                   treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                                   selected.genes = c("CTNNA2", "CCSER1", "FILIP1L", "DHRS3", "PCED1B", #Top 5 downregulated
                                                                      "SLC1A1", "A2ML1-AS1", "CATSPERB", "RIMKLB", "PAEP")) # Top 5 upregulated
    
    SOX9posLGR5neg.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                       celltype.x = "SOX9+ LGR5-", group.name = "Restored",
                                                       treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                                       selected.genes = c("MAP3K5", "A2ML1-AS1", "MALAT1", "GABRP", #Top 5 downregulated
                                                                          "CLDN10", "NEAT1", "RIMKLB", "SLPI", "PAEP")) # Top 5 upregulated
    
    SOX9posLGR5pos.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                       celltype.x = "SOX9+ LGR5+", group.name = "Restored",
                                                       treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                                       selected.genes = c("SEMA3E", "ROBO2", "CCSER1", "CTNNA2", "THUMPD1", #Top 5 downregulated
                                                                          "RIMKLB", "ITGA2", "LAMC2", "ADAMTS9", "PAEP")) # Top 5 upregulated
  } else if (select.idents == "Immune_labelled") {
    
    # Plotting control vs. PCOS
    uNK1.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                celltype.x = "uNK 1", group.name = "Ctrl_PCOS_top10",
                                                treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                                selected.genes = c("OPRM1", "PCBP3", "BACH2", "RIPOR2", "MGAT5", #Top 5 downregulated
                                                                   "ACOXL", "SLPI", "GPX3", "CXCL14", "PAEP")) # Top 5 upregulated
    
    uNK2.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                celltype.x = "uNK 2", group.name = "Ctrl_PCOS_top10",
                                                treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                                selected.genes = c("PCBP3", "DPYD-AS1", "COL3A1", "MGAT5", "MMP7", #Top 5 downregulated
                                                                   "DEFB1", "SLPI", "GPX3", "CXCL14", "PAEP")) # Top 5 upregulated
    
    uNK3.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                celltype.x = "uNK 3", group.name = "Ctrl_PCOS_top10",
                                                treatment.levels = c("Control", "PCOS_W0"),
                                                selected.genes = c("NPAS3", "PAEP")) # Only two DEG's in ctrl vs. PCOS
    
    uM1.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                celltype.x = "uM 1", group.name = "Ctrl_PCOS_top10",
                                                treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                                selected.genes = c("COLEC12", "GPNMB", "KCNMA1", "CYP27A1", "FTL", #Top 5 downregulated
                                                                   "ANPEP", "PID1", "SLC16A10", "NAMPT", "PAEP")) # Top 5 upregulated
    
    uM2.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                celltype.x = "uM 2", group.name = "Ctrl_PCOS_top10",
                                                treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                                selected.genes = c("RBPJ", "AKAP6", "KCNMA1", "GPNMB", "NPL", #Top 5 downregulated
                                                                   "SIPA1L1", "SLPI", "SPP1", "SLC16A10", "PAEP")) # Top 5 upregulated
    
    CD4.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                celltype.x = "T-cells CD4+", group.name = "Ctrl_PCOS",
                                                treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                                selected.genes = c("COL3A1", "RPL17", "TMSB4X", #Top 5 downregulated
                                                                   "SAT1", "IL7R", "BACH2", "CRYBG1", "PDE3B", "FAAH2", # Top 5 upregulated
                                                                   "FOXP3")) # Fill out to get to 10 genes
 
    CD8.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                celltype.x = "T-cells CD8+", group.name = "Ctrl_PCOS_top10",
                                                treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                                selected.genes = c("LGALS1", "FTL", "TMSB4X", "RPL10", #Top 5 downregulated
                                                                   "WFDC2", "NPAS3", # Top 5 upregulated
                                                                   "FOXP3", "IL1B", "KIT", "CD27")) # Fill out to get to 10 genes

    
    
    # Plotting all groups
    uNK1.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                             celltype.x = "uNK 1", group.name = "Restored_Met", scale.min.x = 0, scale.max.x = 60,
                                             treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                             selected.genes = c("OPRM1", "PCBP3", "BACH2", "RIPOR2", "MGAT5", #Top 5 downregulated
                                                                "ACOXL", "SLPI", "GPX3", "CXCL14", "PAEP")) # Top 5 upregulated
    
    uNK2.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                             celltype.x = "uNK 2", group.name = "Restored_Met", scale.min.x = 0, scale.max.x = 60,
                                             treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                             selected.genes = c("DPYD-AS1", "COL3A1", "MMP7", "DTHD1", "OPRM1",  #Top 5 downregulated
                                                                "DEFB1", "SLPI", "GPX3", "CXCL14", "PAEP")) # Top 5 upregulated
    
    uM1.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                            celltype.x = "uM 1", group.name = "Restored_Met", scale.min.x = 0, scale.max.x = 60,
                                            treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                            selected.genes = c("ANKUB1", "DPYD-AS1", "STARD13", "SNHG14", "COL3A1", #Top 5 downregulated
                                                               "SLC43A2", "ANPEP", "SLC16A10", "NAMPT", "PAEP")) # Top 5 upregulated
    
    
    uM2.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                            celltype.x = "uM 2", group.name = "Restored_Met", scale.min.x = 0, scale.max.x = 60,
                                            treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                            selected.genes = c("SIPA1L1", "SPP1", "PAEP")) 
    
    CD4.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                            celltype.x = "T-cells CD4+", group.name = "Restored_Met", scale.min.x = 0, scale.max.x = 60,
                                            treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                            selected.genes = c("RPL17", "IL7R", "CRYBG1", "FAAH2"))
    
    CD8.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                            celltype.x = "T-cells CD8+", group.name = "Restored_Met", scale.min.x = 0, scale.max.x = 60,
                                            treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                            selected.genes = c("WFDC2")) # Top 5 upregulated
    
    # Plotting all groups for uNK1, uNK2, uM1, uM2
    uNK1.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                             celltype.x = "uNK 1", group.name = "Restored_Met_LS", scale.min.x = 0, scale.max.x = 60,
                                             treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                             selected.genes = c("COL3A1", "", "", "", "", #Top 5 downregulated
                                                                "ACOXL", "SLPI", "GPX3", "CXCL14", "PAEP")) # Top 5 upregulated
    
    uNK2.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                             celltype.x = "uNK 2", group.name = "Restored_Met_LS", scale.min.x = 0, scale.max.x = 60,
                                             treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                             selected.genes = c("DPYD-AS1", "COL3A1", "MMP7", "DTHD1", "OPRM1",  #Top 5 downregulated
                                                                "DEFB1", "SLPI", "GPX3", "CXCL14", "PAEP")) # Top 5 upregulated
    
    uM1.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                            celltype.x = "uM 1", group.name = "Restored_Met_LS", scale.min.x = 0, scale.max.x = 60,
                                            treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                            selected.genes = c("ANKUB1", "DPYD-AS1", "STARD13", "SNHG14", "COL3A1", #Top 5 downregulated
                                                               "SLC43A2", "ANPEP", "SLC16A10", "NAMPT", "PAEP")) # Top 5 upregulated
    
    
    uM2.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                            celltype.x = "uM 2", group.name = "Restored_Met_LS", scale.min.x = 0, scale.max.x = 60,
                                            treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"),
                                            selected.genes = c("SIPA1L1", "SPP1", "PAEP")) 
    
  } else if (select.idents == "Endothelial_labelled") {
    Artery.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                             celltype.x = "Endothelial Artery", group.name = "Ctrl_PCOS_top11",
                                             treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                             selected.genes = c("PDE3A", "MGP", "FN1", "SULF1", "CLDN5", "FTL", 
                                                                "FAM118A", "COBLL1", "TSPAN18", "RCAN2", "HDAC9")) # Top 5 upregulated
    
    Vein.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                             celltype.x = "Endothelial Vein", group.name = "Ctrl_PCOS_top10",
                                             treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                             selected.genes = c("MBD5", "TIMP3", "CASP8AP2", "DSN1", "IPO7", "PTMA", "RPS24")) # Top 5 downregulated
    
    
    Proliferative.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                             celltype.x = "Endothelial proliferative", group.name = "Ctrl_PCOS_top10",
                                             treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                             selected.genes = "NPAS3")
    
    Lymphatic.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                            celltype.x = "Lymphatic", group.name = "Ctrl_PCOS_top10",
                                            treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                            selected.genes = c("FN1", "CD9", "DSCAM", "NEAT1", "MGAT4C")) # Top regulated
    
    # Plotting Ctrl vs. PCOS genes in all groups
    Artery.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                               celltype.x = "Endothelial Artery", group.name = "Ctrl_PCOS_top_in_all",
                                               treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                               selected.genes = c("PDE3A", "MGP", "FN1", "SULF1", "CLDN5", "FTL", 
                                                                  "FAM118A", "COBLL1", "TSPAN18", "RCAN2", "HDAC9")) # Top 5 upregulated
    
    Vein.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                             celltype.x = "Endothelial Vein", group.name = "Ctrl_PCOS_top_in_all",
                                             treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                             selected.genes = c("MBD5", "TIMP3", "CASP8AP2", "DSN1", "IPO7", "PTMA", "RPS24")) # Top 5 downregulated
    
    
    Proliferative.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                      celltype.x = "Endothelial proliferative", group.name = "Ctrl_PCOS_top_in_all",
                                                      treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                      selected.genes = "NPAS3")
    
    Lymphatic.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                  celltype.x = "Lymphatic", group.name = "Ctrl_PCOS_top_in_all",
                                                  treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                  selected.genes = c("FN1", "CD9", "DSCAM", "NEAT1", "MGAT4C")) # Top regulated
    
    
    # Plotting genes of interest found in GO of BP
    Artery.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                               celltype.x = "Endothelial Artery", group.name = "Ctrl_PCOS_GO_BP",
                                               treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                               selected.genes = c("PDE3A", "FN1", "SULF1", "CLDN5"))
    
    Vein.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                             celltype.x = "Endothelial Vein", group.name = "Ctrl_PCOS_GO_BP",
                                             treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                             selected.genes = c("MBD5", "TIMP3", "CASP8AP2", "PAEP")) # Top 5 downregulated
    
    Lymphatic.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                  celltype.x = "Lymphatic", group.name = "Ctrl_PCOS_GO_BP",
                                                  treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                  selected.genes = c("FN1", "CD9", "DSCAM", "NEAT1")) # Top regulated
    
  } else if (select.idents == "Stromal_labelled") {
    
    # Plotting Ctrl vs PCOS
    Stroma1.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                               celltype.x = "Stroma 1", group.name = "Ctrl_PCOS_top10",
                                               treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                               selected.genes = c("PIK3C2G", "SFRP1", "SLC8A1-AS1", "SMOC2", "ROBO2", 
                                                                  "BICC1", "MAP3K5", "COL7A1", "CEMIP", "NEAT1")) # Top 5 upregulated
    
    Stroma2.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                             celltype.x = "Stroma 2", group.name = "Ctrl_PCOS_top10",
                                             treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                             selected.genes = c("PIK3C2G", "ROBO2", "SFRP1", "CNTN1", "SLC8A1-AS1", 
                                                                "KCNIP1", "COL7A1", "HPSE2", "MAP3K5", "NEAT1")) # Top 5 downregulated
    
    
    StromaProliferative.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                      celltype.x = "Stroma proliferative", group.name = "Ctrl_PCOS_top10",
                                                      treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                                      selected.genes = c("LRRTM4", "LUM", "KIF26B", "SFRP1", "SULF1", 
                                                                         "NAV2", "KCNIP4", "PTCH1", "HPSE2", "APOLD1"))
    
    Fibroblast.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                  celltype.x = "Fibroblast", group.name = "Ctrl_PCOS_top10",
                                                  treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                                  selected.genes = c("SLC8A1-AS1", "CNTN1", "SLC8A1", "KIF26B", "SMOC2", 
                                                                     "ITGA2", "STC1", "CEMIP", "MMP10", "MMP3")) # Top regulated
    
    uSMC.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated.Ctrl_PCOS,
                                                   celltype.x = "uSMC", group.name = "Ctrl_PCOS_top10",
                                                   treatment.levels = c("Control", "PCOS_W0"), scale.min.x = 0, scale.max.x = 100,
                                                   selected.genes = genes <- c("ROBO2", "SFRP1", "GNA14-AS1", "LIPI", "MMP11", 
                                                                               "MIR4435-2HG", "RUNX1", "PALM2-AKAP2", "CYTOR", "NEAT1")) # Top regulated
    
    # Plotting all groups with Ctrl vs PCOS genes
    Stroma1.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                celltype.x = "Stroma 1", group.name = "Ctrl_PCOS_top10_in_all",
                                                treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                selected.genes = c("PIK3C2G", "SFRP1", "SLC8A1-AS1", "SMOC2", "ROBO2", 
                                                                   "BICC1", "MAP3K5", "COL7A1", "CEMIP", "NEAT1")) # Top 5 upregulated
    
    Stroma2.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                celltype.x = "Stroma 2", group.name = "Ctrl_PCOS_top10_in_all",
                                                treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                selected.genes = c("PIK3C2G", "ROBO2", "SFRP1", "CNTN1", "SLC8A1-AS1", 
                                                                   "KCNIP1", "COL7A1", "HPSE2", "MAP3K5", "NEAT1")) # Top 5 downregulated
    
    
    StromaProliferative.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                            celltype.x = "Stroma proliferative", group.name = "Ctrl_PCOS_top10_in_all",
                                                            treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                            selected.genes = c("LRRTM4", "LUM", "KIF26B", "SFRP1", "SULF1", 
                                                                               "NAV2", "KCNIP4", "PTCH1", "HPSE2", "APOLD1"))
    
    Fibroblast.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                   celltype.x = "Fibroblast", group.name = "Ctrl_PCOS_top10_in_all",
                                                   treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                   selected.genes = c("SLC8A1-AS1", "CNTN1", "SLC8A1", "KIF26B", "SMOC2", 
                                                                      "ITGA2", "STC1", "CEMIP", "MMP10", "MMP3")) # Top regulated
    
    uSMC.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                             celltype.x = "uSMC", group.name = "Ctrl_PCOS_top10_in_all",
                                             treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                             selected.genes = genes <- c("ROBO2", "SFRP1", "GNA14-AS1", "LIPI", "MMP11", 
                                                                         "MIR4435-2HG", "RUNX1", "PALM2-AKAP2", "CYTOR", "NEAT1")) # Top regulated
    
    
    # Plotting restored after PCOS metformin
    Stroma1.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                celltype.x = "Stroma 1", group.name = "PCOS_Met_restored",
                                                treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                selected.genes = c("ROBO2", 
                                                                   "SAT1", "CBLB", "MAML3", "SNX9", "BICC1", "MAP3K5", "COL7A1", "CEMIP", "NEAT1")) # Top 5 upregulated
    
    Stroma2.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                celltype.x = "Stroma 2", group.name = "PCOS_Met_restored",
                                                treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                selected.genes = c("ROBO2", 
                                                                   "COL4A1", "PDE10A", "BICC1", "CBLB", "SNX9", "MAML3", "COL7A1", "MAP3K5", "NEAT1")) # Top 5 downregulated
    
    
    StromaProliferative.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                            celltype.x = "Stroma proliferative", group.name = "PCOS_Met_restored",
                                                            treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                            selected.genes = c("ROBO2", "CCSER1", "CGNL1", "AL049828.1", "NEAT1"))
    
    Fibroblast.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                                   celltype.x = "Fibroblast", group.name = "PCOS_Met_restored",
                                                   treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                                   selected.genes = c("LINC02456", "ROBO2", 
                                                                      "PAEP", "MAP3K5", "INHBA", "NEAT1", "ITGA2", "STC1", "CEMIP", "MMP10", "MMP3")) # Top regulated
    
    uSMC.overlap.plot = DEG.overlap.plotting(seurat.x = endo.integrated,
                                             celltype.x = "uSMC", group.name = "PCOS_Met_restored",
                                             treatment.levels = c("Control", "PCOS_W0", "PCOS_W16_Met", "PCOS_W16_LS"), scale.min.x = 0, scale.max.x = 100,
                                             selected.genes = c("DLEU1", 
                                                                "PRR16", "COL7A1", "MAP3K5", "PTPRE", "GLIS3", "CBLB", "LOXL2", "RUNX1", "NEAT1")) # Top regulated
    
    
  }
  
} else if (server.run == FALSE) {
  
  #input.dir.list = list.files(path = paste0(Input.dir, celltype), pattern = "*_sign_DEG.xlsx", full.names = TRUE)

  celltype.list = list.dirs(Input.dir, full.names = FALSE)
  celltype.list = celltype.list[nchar(celltype.list) != 0]
  celltype.list = celltype.list[!str_detect(celltype.list,pattern="/GSEA_enrichment_analysis")]
  
  Overlap.output = lapply(celltype.list, Overlap.DEGs)

  # Create a woerkbook to add worksheets of GO results to
  Overlap.output = lapply(celltype.list[6], GO_enrichment_analysis, skip.overlap = skip.overlap.flag)

}
