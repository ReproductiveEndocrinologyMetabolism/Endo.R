    _______ .__   __.  _______   ______      .______      
    |   ____||  \ |  | |       \ /  __  \     |   _  \     
    |  |__   |   \|  | |  .--.  |  |  |  |    |  |_)  |    
    |   __|  |  . `  | |  |  |  |  |  |  |    |      /     
    |  |____ |  |\   | |  '--'  |  `--'  |  __|  |\  \----.
    |_______||__| \__| |_______/ \______/  (__) _| `._____|
                                               

README
Author: Gustaw Eriksson  
Date: 01-11-2024  
Version: 1.0   
Contact: gustaw.eriksson@ki.se  

# Introduction
Endometrium analysis in R (Endo.R) is a downstream data analysis pipeline based on Seurat developed to analysis CellRanger output for used to study 10x Genomics based single-nuclei (sn) RNA-seq on endometrium tissue. The pipeline has been developed to handle datasets of large numbers of cells from multiple study groups, specifically control vs. case vs. treatment. Depending on the size of the dataset, high performance cluster environment has to be used to run the pipeline. Endo.R was used in the following project and publication:  

Eriksson, G., Li, C., et al. Single-Cell Profiling of the Human Endometrium in Polycystic Ovary Syndrome: Uncovering Disease Signatures and Treatment Responses. Manuscript (2024)

The analysed 10X snRNA-seq data is archived in the European Genome-phenome Archive (EGA), ID: EGAD50000001017 and is released by request.

For validation of the snRNA-seq data. spatial RNA-sequencing (Stereo-seq) was performed. This repository also contains code used for the analysis. 

## Description
Below are instructions on how to build the project directory required for the snRNA-seq and Stereo-seq pipeline, required input and description about the scripts.

# 10X snRNA-seq of endometrium
## Project directory
Create a project directory and copy the scripts to it.

## Input directory
Create a data directory in the project directory as such: EndoR/Data/. Within the /Data directory, create a directory named /CellRanger_Count.  

In Data/CellRanger_Count, create a folder for each sample (Data/CellRanger_Count/Sample_X) containing Cell Ranger output for each sample and the required "filtered_feature_bc_matrix" folder which is generated after running Cell Ranger.    Create the directory by running:

```
mkdir EndoR  
mkdir EndoR/Data  
mkdir EndoR/Data/CellRanger_Count  
cp -r CellRanger_sample_X EndoR/Data/CellRanger_Count
```
Or by running Create_directory.sh script in the project directory:

```
bash Create_directory.sh
```
## Output directory
The scripts in the pipeline will automatically generate a /Output directory in the project directory. Within it, output folders for each script and analysis will be outputted.

## Running the scripts
Now you are ready to run the scripts. Run the scripts in the numbered order and make sure to manually check the output between each scripts. Below follows brief description of each script and step.

### 0. Individual sample preprocessing
Each Cell Ranger (filtered feature bc matrix) is loaded to generate a Seurat object. The object is then filtered according to set QC cutoffs for mitochondrian, ribosomal and hemoglobin %RNA as well as number of feature, features per cell to remove any nuclei of low quality. All filtration step output plots and data in the generated output folder. The object can thereafter be SCT transformed for sample integration.

### 1. Sample integration
All preprocessed sample are integrated to one uniform object. In this script, SCT transformation is applied but the script can be edited to use other integration methods. Moreover, the script saves the integration anchors and outputs plots related to the QC.

To be noted, the script requires manual editing, assigning each sample (orig.ident) to a specific group ("Group_Stage). 

### 2. Main cell type clustering
The integrated Seurat object is clustered and UMAPs with gene markers are generated. Please edit the script to add and remove gene markers to fit your analysis. If requested, automated annotation using SingleR can be applied.

### 3. Main cluster annotation
Based on findings in step 2, the Seurat clusters are annotated to a specific main cell type. In this project, clusters were annoteded as either epithelium, stromal, smooth muscle cells, lymphoid, myeloid, endothelial or lymphatic. Another set of QC plots are generated, grouped based on the main cell type annotation. A new Seurat object is saved including the annotation data.

### 4. Subsetting of the main clusters
The annotated Seurat object is subsetted based on main cell type annotation, generating four new Seurat objects; Epithelium, stroma + smooth muscle cells, immune, and endothelial + lymphatics.  

### 5. Sub-clustering and annotation
Each script is used for one subset from step 4. The Seurat object is split based on sample and re-integrated and re-clustered. Therafter, gene markers for subclusters of the main cell type is projected to identify subclusters which are then annotated. The annoated subset Seurat object is saved and new QC plots are generated.

### 6. Calculating differential gene expression 
Differential gene expression (DGE) analysis can be performed with either Seurat default Wilcoxon Rank Sum Test or MAST hurdle model used in the project to account for sample variance. DGE is performed between control vs. PCOS, PCOS vs. Metformin and PCOS vs. Lifestyle in all subcluster of each Seurat subset. The script generated plots of the differentially expressed genes (DEGs) and tables of the DEGs. 

### 7. Analysis of differential gene expression and gene ontology enrichment
Gene Ontology Enrichment analysis with multiple methods and reference on generated DEG tables. In 07.1, the GOEA is performed and tables of the results are generated. The script generated a table of DEGs of all comparisons in step 6, including information if the DEG is up- or downexpressed, and if it is restored or not in the treatment groups. GOEA is performed on three sets of DEGs; All, upexpressed and downexpressed. Please change p- and q-values if needed. 

The tables are generated in .xlsx allowing for curation of the tables, preparing the data for additional visualisation in the following 07.2 script. Here, GO enrichment is plotted using barpltos visualising if it is restored or not in treatment. 

### 8. Visualisation of genes and ontologies
Genes of interests identified in step 6 and 7 can be further visualised in this step using heatmaps and violinplots. To reduce memory load, the objects can be averaged and saved as new objects of decreased memory. Averaged objects can be used for heatmap generation. 

### 9. Subcluster projection on main Seurat object
To perform step 10-13, the original main seurat object of all integrated samples and clusters, is re-annoted with the identified subclusters of all main cell ttypes. This is performed by combining the annotations per UMI of all subclusters to one combined vector and adding it to the main seurat object. The re-annotated and re-combined Seurat object is saved and is used in the below steps. 

### 10. CellChat ligand-receptor interaction analysis
First in 10.1, the re-combined Seurat object from step 9 is loaded and subset into the four group of which a CellChat object is generated by each. The CellChat objects are saved and combinded into one used in 10.2 when comparing the groups.

In 10.2, interactions strengts of all subclusters are compared and visualised to determine if any differences can be identified between the groups. In the following script 10.3, all subclusters are compared for interacting ligand-receptors and associated pathways which are presented in generated tables. 

The DEG and CellChat tables are combined and overlapped in 10.4, filtering the CellChat tables for potential DEG which then are visualised with Dotplots of ligands and receptors of the compared groups. The script allows for selecting which pathways and DEGs to be visualised by curating the CellChat table. 

### 11. CELLECT - CELLEX GWAS integration
TINA, PLEASE ADD

### 12. Correlation analysis of gene expression and clinical measurments


### 13. Cluster proportion analysis

### 14. Folder of additional scripts for visualisation

# Stereo-seq of endometrium
## Project directory
Create a project directory and copy the scripts to it.

## Input directory
Create a data directory in the project directory as such: EndoR/Data/. Within the /Data directory, create a directory named /CellRanger_Count.  

In Data/CellRanger_Count, create a folder for each sample (Data/CellRanger_Count/Sample_X) containing Cell Ranger output for each sample and the required "filtered_feature_bc_matrix" folder which is generated after running Cell Ranger.    Create the directory by running:

```
mkdir EndoR  
mkdir EndoR/Data  
mkdir EndoR/Data/CellRanger_Count  
cp -r CellRanger_sample_X EndoR/Data/CellRanger_Count
```
Or by running Create_directory.sh script in the project directory:

```
bash Create_directory.sh
```
## Output directory
The scripts in the pipeline will automatically generate a /Output directory in the project directory. Within it, output folders for each script and analysis will be outputted.

## Running the scripts
Now you are ready to run the scripts. Run the scripts in the numbered order and make sure to manually check the output between each scripts. Below follows brief description of each script and step.

### 0. Individual sample preprocessing

### 1. UMAP label transfer annotation

### 2. Annotation quality check and filtering

### 3. Sub-clustering annotated spatial object

### 4. Marker visualisation on full and targeted spatial object

# Acknowledgement

The computations were enabled by resources provided by the Swedish National Infrastructure for Computing (SNIC) at UPPMAX partially funded by the Swedish Research Council through grant agreement no. 2018-05973.
