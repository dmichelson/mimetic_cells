# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 09:46:10 2020

@author: dmich
"""
#%%
import os
import pandas as pd
import scanpy as sc
import numpy as np
import scipy as scp
import sklearn
import matplotlib.pyplot as plt
import matplotlib
import sys
import loompy
import scipy.optimize
import scvelo as scv
import glob
import pickle
import anndata
from collections import Counter
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform

scv.settings.set_figure_params('scvelo')
plt.rcParams['figure.dpi'] = 600

#%%
## read loom
os.chdir("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/velocity")

adata = anndata.read_loom("./loom/post_meclo_adata_2020-11-09.loom")

# plt.rcParams['image.cmap'] = 'tab20'
# tab14 = ['#59A14F','#FF9D9A','#499894','#FFBE7D','#E15759','#F1CE63','#BAB0AC','#79706E',
#          '#86BCB6','#F28E2B','#8CD17D', '#B6992D','#4E79A7', '#A0CBE8']
# adata.uns['Clusters_colors'] = tab14

## scvelo preprocess
# umap = pd.read_csv("/n/groups/cbdm-db/dam41/postaire_meclo_scrna/seurat-out/cell_embeddings.csv")
# umap = umap.iloc[:,1:]
# adata.obsm['X_umap'] = umap.values

#%%
os.chdir("C:/Users/dmich/Dropbox (MIT)/CBDM Lab MIT/Data/MEC_scRNA/DAM_082620/velocity")
adata = anndata.read_loom("./loom/post_meclo_adata_2020-11-09.loom")
transcripts = pd.read_csv("./text_inputs/gene_names.csv")
adata.var.index = transcripts.values[:,1]

#%%
cell_clusters = pd.read_csv("./text_inputs/active_clusters.csv")
cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'Cell ID'})
adata.obs["Clusters"] = cell_clusters["x"].values

idx = ( adata.obs['Clusters']=="Gut/Liver") | ( adata.obs['Clusters']=="Mcell" )
# idx = ( adata.obs['Clusters']=="Neuroendocrine" ) | ( adata.obs['Clusters']=="Skin" ) | (adata.obs['Clusters']=="Aire-stage")
# idx = ( adata.obs['Clusters']=="Neuroendocrine" ) | ( adata.obs['Clusters']=="Ciliated" ) | ( adata.obs['Clusters']=="Ciliated" )
# idx = ( adata.obs['Clusters']=="Muscle" ) | ( adata.obs['Clusters']=="Tuft2" ) | ( adata.obs['Clusters']=="Tuft1" )| ( adata.obs['Clusters']=="Skin, basal" )

subset_adata = adata[idx].copy()
subset_adata.var_names_make_unique()

#%%
# sc.pp.calculate_qc_metrics(subset_adata)
# sc.pp.filter_genes(subset_adata, min_cells=3)
# sc.pp.log1p(subset_adata)
# sc.pp.highly_variable_genes(subset_adata)
# sc.pp.normalize_total(subset_adata)
# sc.pp.pca(subset_adata, n_comps=50)
# sc.pp.neighbors(subset_adata, n_neighbors=20, n_pcs=5)

scv.pp.filter_and_normalize(subset_adata, min_shared_counts=20, n_top_genes=20000)
scv.pp.moments(subset_adata, n_pcs=5, n_neighbors=30)


#%%
sc.tl.umap(subset_adata)
sc.tl.diffmap(subset_adata)
# sc.tl.leiden(subset_adata, resolution=1.5)
# sc.tl.embedding_density(subset_adata, basis='umap')

sc.pl.umap(subset_adata, color='Clusters')
sc.pl.diffmap(subset_adata, color='Clusters')
# sc.pl.embedding_density(subset_adata, basis='umap', save="_subset_density_2021-04-09.pdf")

scv.tl.recover_dynamics(subset_adata, n_jobs=6)

scv.tl.velocity(subset_adata, mode="dynamical")
scv.tl.velocity_graph(subset_adata)
scv.pl.velocity_embedding(subset_adata, basis="diffmap", color="Clusters", add_outline=['Gut/Liver','Mcell'],
                          outline_width=[0.1,0.01], arrow_length=2, arrow_size=2, linecolor="black", legend_loc="upper right"), 
                          save="gut_mcell_arrow_2021-06-16.pdf")
scv.pl.velocity_embedding_grid(subset_adata, basis="diffmap", color="Clusters", add_outline=['Gut/Liver','Mcell'],
                          outline_width=[0.1,0.01], size=50, legend_loc="upper right", save="gut_mcell_grid_2021-06-16.pdf")
scv.pl.velocity_embedding_stream(subset_adata, basis="diffmap", color="Clusters", add_outline=['Gut/Liver','Mcell'],
                          outline_width=[0.1,0.01], size=50, legend_loc="upper right"), save="gut_mcell_stream_2021-06-16.svg")


scv.pl.velocity_embedding(subset_adata, basis="umap", color="Clusters", legend_loc='right')
scv.pl.velocity_embedding(subset_adata, basis="diffmap", color="Clusters", legend_loc='right')
scv.pl.velocity_embedding_grid(subset_adata, basis="umap", color="Clusters", legend_loc='right')
scv.pl.velocity_embedding_grid(subset_adata, basis="diffmap", color="Clusters", legend_loc='right')
scv.pl.velocity_embedding_stream(subset_adata, basis="umap", color="Clusters", legend_loc='right')
scv.pl.velocity_embedding_stream(subset_adata, basis="diffmap", color="Clusters", legend_loc='right')

mcell_sig = pd.read_csv("./text_inputs/mcell_sig_2021-06-16.csv")
gut_sig = pd.read_csv("./text_inputs/enterocyte_mature_prox_sig_2021-06-16.csv")
liver_sig = pd.read_csv("./text_inputs/hepatocyte_sig_2021-06-16.csv")

gut_sig.squeeze()

gut_sig = gut_sig.rename(columns={'na.omit(intestinal_sig$Enterocyte.Mature.Proximal)':'V1'})
liver_sig = liver_sig.rename(columns={'na.omit(hepatocyte_sig)':'V1'})

gutliver_sig = gut_sig.append(liver_sig, ignore_index=True)

sc.tl.score_genes(subset_adata, liver_sig.squeeze())
subset_adata.obs.score
sc.pl.diffmap(subset_adata, color=['score'], title="Gut/liver", add_outline=['Flg'], 
              outline_width=[0.1,0.01], vmin="p5", vmax="p95", color_map="Purples" )

sc.pl.diffmap(subset_adata, color=['Hamp'], add_outline=['Gut/Liver','Mcell'],
              outline_width=[0.1,0.01], vmin="p5", vmax="p95", color_map="Purples")

sc.pl.diffmap(subset_adata, color=['Aldob','Saa3','Pigr','Vil1','Marcksl1','Tnfaip2','Gp2','Ccl9','Ccl20'], add_outline=['Gut/Liver','Mcell'],
              outline_width=[0.1,0.01], vmin="p5", vmax="p95", color_map="Purples", ncols=3, save="_gut_mcell_markers_2021-06-16.png")

#%%
cell_clusters = pd.read_csv("./text_inputs/active_clusters.csv")
cell_clusters = cell_clusters.rename(columns={'Unnamed: 0': 'Cell ID'})
adata.obs["Clusters"] = cell_clusters["x"].values

# idx = ( adata.obs['Clusters']=="Gut/Liver") | ( adata.obs['Clusters']=="Mcell" )
idx = ( adata.obs['Clusters']=="Skin, keratinized" ) | ( adata.obs['Clusters']=="Skin, basal" )
# idx = ( adata.obs['Clusters']=="Neuroendocrine" ) | ( adata.obs['Clusters']=="Ciliated" ) | ( adata.obs['Clusters']=="Ciliated" )
# idx = ( adata.obs['Clusters']=="Muscle" ) | ( adata.obs['Clusters']=="Tuft2" ) | ( adata.obs['Clusters']=="Tuft1" )| ( adata.obs['Clusters']=="Skin, basal" )

subset_adata = adata[idx].copy()
subset_adata.var_names_make_unique()

#%%

scv.pp.filter_and_normalize(subset_adata, min_shared_counts=20, n_top_genes=20000)
scv.pp.moments(subset_adata, n_pcs=5, n_neighbors=30)


#%%
sc.tl.umap(subset_adata)
sc.tl.diffmap(subset_adata)
sc.tl.leiden(subset_adata, resolution=0.5)
sc.tl.rank_genes_groups(subset_adata, groupby="leiden")
# sc.tl.embedding_density(subset_adata, basis='umap')

sc.pl.umap(subset_adata, color='leiden')
sc.pl.diffmap(subset_adata, color='leiden')
# sc.pl.embedding_density(subset_adata, basis='umap', save="_subset_density_2021-04-09.pdf")

scv.tl.recover_dynamics(subset_adata, n_jobs=6)
scv.tl.velocity(subset_adata, mode="dynamical")
scv.tl.velocity_graph(subset_adata)

scv.pl.velocity_embedding(subset_adata, basis="diffmap", color="black",
                          outline_width=[0.1,0.01], arrow_length=2, arrow_size=2, linecolor="black", legend_loc="upper right",
                          save="skin_arrow_2021-06-16.png")

scv.pl.velocity_embedding_grid(subset_adata, basis="diffmap", color="black",
                          outline_width=[0.1,0.01], size=50, legend_loc="upper right" ,
                          save="skin_grid_2021-06-16.png")

scv.pl.velocity_embedding_stream(subset_adata, basis="diffmap", color="black",
                          outline_width=[0.1,0.01], size=50, legend_loc="upper right" ,
                          save="skin_stream_2021-06-16.png")


sc.pl.diffmap(subset_adata, color=['Krt10'], add_outline=['Flg'],
              outline_width=[0.1,0.01], vmin="p5", vmax="p95", color_map="Purples")

sc.pl.diffmap(subset_adata, color=['Krt5','H2-Aa','Itgb1','Ivl','Flg','Il1f5'], add_outline=['Skin, basal','Skin, keratinized'],
              outline_width=[0.1,0.01], vmin="p5", vmax="p95", color_map="Purples", ncols=3,
              save="_skin_markers_2021-06-16.png")

skin_basal = pd.read_csv("./text_inputs/skin_basal_2021-06-16.txt")
skin_diff = pd.read_csv("./text_inputs/skin_differentiated_2021-06-16.txt")
skin_keratinized = pd.read_csv("./text_inputs/skin_keratinized_2021-06-16.txt")

sc.tl.score_genes(subset_adata, skin_keratinized.squeeze())
subset_adata.obs.score
sc.pl.diffmap(subset_adata, color=['score'], title="Keratinized", add_outline=['Flg'], 
              outline_width=[0.1,0.01], vmin="p5", vmax="p95", color_map="Purples",
              save="_skin_keratin_sig_2021-06-16.png")


