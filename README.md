                                       _______ .__   __.  _______   ______      .______      
                                      |   ____||  \ |  | |       \ /  __  \     |   _  \     
                                      |  |__   |   \|  | |  .--.  |  |  |  |    |  |_)  |    
                                      |   __|  |  . `  | |  |  |  |  |  |  |    |      /     
                                      |  |____ |  |\   | |  '--'  |  `--'  |  __|  |\  \----.
                                      |_______||__| \__| |_______/ \______/  (__) _| `._____|
                                                                                     

README
Author: Gustaw Eriksson  
Date: 01-11-2024  
Version: 1.0. 
Contact: gustaw.eriksson@ki.se  

# Introduction
Endometrium analysis in R (Endo.R) is a downstream data analysis pipeline based on Seurat developed to analysis CellRanger output for used to study 10x Genomics based single-nuclei (sn) RNA-seq on endometrium tissue. The pipeline has been developed to handle datasets of large numbers of cells from multiple study groups, specifically control vs. case vs. treatment. Depending on the size of the dataset, high performance cluster environment has to be used to run the pipeline. Endo.R was used in the following project and publication:  
Eriksson, G., Li, C., et al. Single-Cell Profiling of the Human Endometrium in Polycystic Ovary Syndrome: Uncovering Disease Signatures and Treatment Responses. Manuscript (2024)

# Description
Below are instructions on how to build the project directory required for the pipeline, required input and description about the scripts

## Project directory
Create a project directory and copy the scripts to it.

## Input directory
Create a data directory in the project directory as such: MUSENDIPOSER/Data/. Within the /Data directory, create a directory named /CellRanger_Count.  

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
Now you are ready to run the scripts. Run the scripts in the numbered order and make sure to manually check the output between each scripts.

# Acknowledgement

The computations were enabled by resources provided by the Swedish National Infrastructure for Computing (SNIC) at UPPMAX partially funded by the Swedish Research Council through grant agreement no. 2018-05973.
