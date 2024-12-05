#!/usr/bin/env python
# coding: utf-8

# Code to preprocess Stereo .gef data to specifically set the bin size before saving it as an .h5ad object for further downstream analysis in R and Seurat. Before starting, set the kernel to the "st" kernel allowing to run Stereopy package. The script is also ment to be runned on the in-house Deng-server or similar HPC environment.
# 
# The script and pre-processing analysis is performed on one sample at the time. Therefore, the path to the specific .gef file has to be manually changed below.

# In[72]:


#import stereo as st
#import warnings

import sys
import stereo as st
import pandas as pd
import numpy as np
import time
import warnings
warnings.filterwarnings('ignore')


# Sanity check to see if Stereopy is loaded and what version it is

# In[73]:


st.__version__


# Vector are generated including key marker genes and genes of interest

# In[74]:


epithelium_subtype = ['PTGS1', 'VTCN1', 'SLC26A7', 'LGR5', 'KRT5', 'WNT7A', 'CPM', 'IHH', 'EMID1', 'PPARG', 'C2CD4A', 'SLC18A2',
                     'PAEP', 'CXCL14', 'MKI67', 'HMGB2', 'AR', 'CDC20B', 'CCNO', 'HES6']
stroma_subtype = ['ESR1', 'PGR', 'IGF1', 'ECM1', 'PAEP', 'OGN', 'TOP2A', 'MKI67', 'THY1', 'COL1A1', 'PCOLCE', 'C7', 'ACTA2', 'ACTG2', 'MCAM']
all_marker_genes = ['EPCAM', 'CPM', 'LGR5', 
                    'IGF1', 'DCN', 'COL6A1', 
                   'GUCY1A2', 'ACTA2', 'NOTCH3',
                   'CD14', 'CSF1R', 'LYZ',
                   'STK17B', 'NCAM1', 'CCL5', 'CD2',
                   'PCDH17', 'VWF',
                   'PROX1', 'FLT4']
print(epithelium_subtype)
print(stroma_subtype)
print(all_marker_genes)


# ## Loading the data
# Setting the path below with above paths

# In[75]:


data_path = './Data/Endometrium_StereoSeq_Spatial_preprocessed_data_230831/C02132D6_027_0_LS/01.StandardWorkflow_Result/GeneExpMatrix/C02132D6.tissue.gef'
st.io.read_gef_info(data_path)
output_bin = './Output/0_Preprocessing_unfiltered_bin30/C02132D6_8_PCOS_unfiltered_bin30.h5ad'
output_log2 = './Output/0_Preprocessing_bin30_filtered_log2/C02132D6_8_PCOS_bin30_log2.h5ad'
output_SCT = './Output/0_Preprocessing_bin30_filtered_SCT/C02132D6_08_PCOS_bin30_SCT.h5ad'
output_ann = './Output/0_Preprocessing_bin30_filtered_SCT_Ann/C02132D6_8_PCOS_bin30_SCT_Ann.h5ad'


# Load the data and generate a StereoExpData object

# In[76]:


data = st.io.read_gef(file_path=data_path, bin_size=30)

# simply type the varibale to get related information
data


# ## Quality Control of non-normalised data
# Generate basic QC measures such as total counts per bin, genes by count and mitochondria percentage. Stomics describes it as:  
# 
# total_counts - the total counts per cell;
# 
# n_genes_by_counts - the number of genes expressed in count maxtrix;
# 
# pct_countss_mt - the percentage of counts in mitochondrial genes.

# In[77]:


data.tl.cal_qc()
data.plt.violin()


# Project the QC data on the spatial scatter

# In[78]:


data.plt.spatial_scatter()


# ## Filtering the non-normalised data
# Stereopy works with three methods to filter the data:  
# data.tl.filter_cells, data.tl.filter_genes, data.tl.filter_coordinates.
# 
# The data can therefore be filtered based on three different levels; cell, gene and coordinate where a bin unit could be treated as a made-up cell provisionally.
# 
# We filter bin units (bin_size is set to 50 at the beginning of our example) based on quality control indicators which have been calculated in QC part. Beforehand, observe the distribution of cells according to scatter plots.

# In[79]:


data.plt.genes_count()


# We filter away bins with >20% mitochondria and with few counts and genes expressed.

# In[80]:


data.tl.filter_cells(
        min_gene=20,
        min_n_genes_by_counts=250,
        pct_counts_mt=20,
        inplace=True
        )
data
data.tl.cal_qc()
data.plt.violin()
data.plt.genes_count()
data.plt.spatial_scatter()


# Checkpoint to save the raw data to be able to redo the filtering

# In[81]:


data.tl.raw_checkpoint()
data.tl.raw


# To retrieve the raw data, run data.tl.reset_raw_data().

# ### Saving the binned data
# Before normalisation, we save the binned and filtered data.

# In[82]:


info = help(st.io.write_h5ad)
print(info)


# In[83]:


st.io.write_h5ad(data, output=output_bin)


# ## Normalization with Log2, currently not used
# We run log2 normalisation to have it in the object. Labelling and annotation is however done with scTransformed counts

# In[84]:


info = help(data.tl.raw_checkpoint)

print(info)


# In[85]:


info = help(data.tl.reset_raw_data)

print(info)


# In[86]:


info = help(data.tl.normalize_total)
print(info)


# In[87]:


info = help(data.tl.log1p)
print(info)


# In[88]:


##data.sparse2array()

##gmean = np.exp(np.log(data.exp_matrix.T + 1).mean(1)) - 1

# preprocessing
##data.tl.raw_checkpoint()
##data.tl.normalize_total(target_sum=1e4)
##data.tl.log1p()

##log_normalize_result = pd.DataFrame([gmean, data.exp_matrix.T.var(1)], index=['gmean', 'log_normalize_variance'], columns=data.gene_names).T

##from stereo.algorithm.sctransform.plotting import plot_log_normalize_var

##fig1=plot_log_normalize_var(log_normalize_result)


# Save the log2 normalised object

# In[89]:


#st.io.write_h5ad(data, output=output_log2)


# ## Normalization via scTransform
# We proceed with normalizing and clustering with scTransform. The default scTransform settings are used, sampling 5000 cells and 2000 genes for estimating parameters. 
# 
# First we perform preprocessing

# In[90]:


info = help(data.tl.sctransform)
print(info)


# ### Restore the data to raw before scTransformation (after log2 normalisation, however this is not working)

# In[91]:


#data.tl.reset_raw_data
#data.tl.raw


# In[92]:


# Preprocessing
#data.tl.sctransform(res_key='sctransform', inplace=False, filter_hvgs=True, n_cells=5000, n_genes=2000)
#data.tl.sctransform(res_key='sctransform', inplace=True, filter_hvgs=False, n_cells=5000, n_genes=2000)
data.tl.sctransform(res_key='sctransform', inplace=False, filter_hvgs=True, n_cells=5000, n_genes=2000)

from stereo.algorithm.sctransform.plotting import plot_residual_var

fig2=plot_residual_var(data.tl.result['sctransform'])


# We embed the data and run PCA and UMAP

# In[93]:


# Embedding
data.tl.pca(use_highly_genes=False, hvg_res_key='highly_variable_genes', n_pcs=20, res_key='pca', svd_solver='arpack')
data.tl.neighbors(pca_res_key='pca', n_pcs=30, res_key='neighbors', n_jobs=8)
data.tl.umap(pca_res_key='pca', neighbors_res_key='neighbors', res_key='umap', init_pos='spectral', spread=2.0)


# In[94]:


# Clustering
data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden')


# Generating a UMAP of the SCT clustering

# In[95]:


data.plt.umap(res_key='umap', cluster_key='leiden')


# We extract the top 50 feature genes

# In[96]:


sct_high_genes = data.tl.result['sctransform'][1]['top_features']
sct_high_genes.tolist()[:50]


# ## Gene distribution
# Investigate how marker genes distribute in the UMAP and spatial dimension.
# 
# First, the UMAP distribution

# In[97]:


#data.plt.umap(res_key='umap', gene_names=epithelium_genes)
#data.plt.umap(res_key='umap', gene_names=stroma_genes)
data.plt.umap(res_key='umap', gene_names=all_marker_genes)


# Checking for the genes in the spatial distribution

# In[98]:


info = help(data.plt.spatial_scatter_by_gene)
print(info)

data.plt.spatial_scatter_by_gene(gene_name=all_marker_genes)


# In[99]:


data.plt.spatial_scatter_by_gene(gene_name=epithelium_subtype)


# In[100]:


data.plt.spatial_scatter_by_gene(gene_name=stroma_subtype)


# ## Identify highly variable genes
# We identify highly variable genes in the bins with their parameters and plot:

# In[101]:


data.tl.highly_variable_genes(
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_top_genes=2000,
        res_key='highly_variable_genes'
        )
data.plt.highly_variable_genes(res_key='highly_variable_genes')


# The gene count is scaled.

# In[102]:


data.tl.scale(max_value=10, zero_center=True)


# ## Spatial hotspots

# In[103]:


info = help(data.tl.spatial_hotspot)
print(info)
data.tl.spatial_hotspot(
                    use_highly_genes=True,
                    use_raw=True,
                    hvg_res_key='highly_variable_genes',
                    model='normal',
                    n_neighbors=30,
                    n_jobs=20,
                    fdr_threshold=0.05,
                    min_gene_threshold=10,
                    res_key='spatial_hotspot',
                    )


# In[104]:


data.plt.hotspot_local_correlations()


# In[105]:


data.plt.hotspot_modules()


# Extract the results of the modules

# In[106]:


data.tl.result['spatial_hotspot']
data.tl.result['spatial_hotspot'].modules


# ## Leiden clustering
# STOmics and Stereopy default output is in Leiden and is therefore used for clustering:

# In[107]:


data.plt.cluster_scatter(res_key='leiden')


# We can perform partial clustering on the data, focusing on certain groups.

# In[108]:


data.plt.umap(res_key='umap', cluster_key='leiden')


# In Jupyter, we can interact with the plot:

# In[109]:


data.plt.interact_cluster(res_key='leiden')


# We can also display the spatial neighbors

# In[110]:


#data.tl.leiden(neighbors_res_key='spatial_neighbors', res_key='spatial_leiden')


# In[111]:


#data.plt.cluster_scatter(res_key='spatial_leiden')


# ## Louvain clustering
# We test Louvain clustering. To reduce the object size, we do not perform it.

# In[112]:


data.tl.louvain(neighbors_res_key='neighbors', res_key='louvain')


# In[113]:


data.plt.cluster_scatter(res_key='louvain')


# In[114]:


data.plt.umap(res_key='umap', cluster_key='louvain')


# ## Phenograph clustering
# We test phenograph clustering. To reduce the object size, we do not perform it.

# In[115]:


#data.tl.phenograph(phenograph_k=30, pca_res_key='pca', res_key='phenograph')


# In[116]:


#data.plt.cluster_scatter(res_key='phenograph')


# ## Find Marker Genes
# Using the Leiden clustering, we can check for marker genes as in Seurat. Default setting is t_test but we change to wilcoxon_test to match Seurat.

# In[117]:


info = help(data.tl.find_marker_genes)
print(info)

data.tl.find_marker_genes(
        cluster_res_key='leiden',
        method='wilcoxon_test',
        use_highly_genes=False,
        use_raw=True
        )


# Plot the marker genes found above.

# In[118]:


data.plt.marker_genes_text(
        res_key='marker_genes',
        markers_num=10,
        sort_key='scores'
        )


# We can generate a dotplot of the data with the top 5 marker genes

# In[119]:


info = help(data.plt.marker_genes_scatter)
print(info)
data.plt.marker_genes_scatter(res_key='marker_genes', markers_num=5)


# ## Save the data before annotation
# Saving the data before annotation with the flavor set to Seurat, allowing for analysis in R with Seurat. 

# In[120]:


#adata = st.io.stereo_to_anndata(data, flavor='seurat', output='./Output/0_Preprocessing_bin30_filtered_SCT/Filename_bin30_SCT.h5ad')
st.io.write_h5ad(data, output=output_SCT)


# ## Cell type annotation
# With the marker genes, we can annotate the object.

# In[121]:


info = help(data.plt.marker_genes_scatter)
print(info)


# In[122]:


data.plt.marker_genes_scatter(res_key='marker_genes', genes=all_marker_genes)
data.plt.marker_genes_scatter(res_key='marker_genes', genes=epithelium_subtype)
data.plt.marker_genes_scatter(res_key='marker_genes', genes=stroma_subtype)
data.plt.umap(res_key='umap', cluster_key='leiden')


# We can add the annotation data in a seperate vector:

# In[123]:


annotation_dict = {
    '1':'Undefined', '2':'Epithelium',
    '3':'Epithelium', '4':'Stroma',
    '5':'Stroma', '6':'Stroma',
    '7':'Undefined', '8':'Epithelium',
    '9':'Stroma', '10':'Epithelium'}


# In[124]:


data.plt.interact_annotation_cluster(
            res_cluster_key='leiden',
            res_marker_gene_key='marker_genes',
            res_key='leiden_annotation'
            )


# In[125]:


data.tl.annotation(
        annotation_information=annotation_dict,
        cluster_res_key='leiden',
        res_key='anno_leiden'
        )


# Plot the annotated data:

# In[126]:


data.plt.cluster_scatter(res_key='anno_leiden')


# ## Saving to .h5ad
# Transform the data to AnnData and save it as an .H5ad file which can later be transformed to .H5seurat for analysis in R in harmony with the 10x data.

# In[127]:


adata = st.io.stereo_to_anndata(data, flavor = 'seurat', output=output_ann)


# In[ ]:




