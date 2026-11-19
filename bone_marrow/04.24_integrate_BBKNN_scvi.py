print('1. scVI-latent BBKNN Setup')
# BBKNN on the scVI latent (04.18, no species key) instead of PCA
import sys
import sklearn.metrics
import sklearn.neighbors
sklearn.neighbors.DistanceMetric = sklearn.metrics.DistanceMetric
sys.modules['sklearn.neighbors.DistanceMetric'] = sklearn.metrics.DistanceMetric
import scanpy as sc
import bbknn
import pandas as pd
import numpy as np
import os
from datasets import selected

for name, d in selected():
    print(f'Load data - {name}:')
    outpath = d['out']
    adata = sc.read_h5ad(d['concat'])
    Z = pd.read_csv(f'{outpath}scvi_latent.csv', index_col = 0)
    adata.obsm['X_scvi_latent'] = Z.loc[adata.obs_names].to_numpy(dtype=float)

    print('Preprocessing')
    adata.obsm['X_raw_pca'] = adata.obsm['X_pca'].copy()
    adata.obsm['X_raw_umap'] = adata.obsm['X_umap'].copy()

    print('2. Run scVI-latent BBKNN')
    bbknn.bbknn(adata, batch_key='orig.ident', use_rep = 'X_scvi_latent')

    print('3. BBKNN Postprocessing')
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_bbknn_scvi'] = adata.obsm['X_umap']

    print('4. scVI-latent BBKNN Save data')
    for df in (adata.obs, adata.var):
            for col in df.columns:
                if isinstance(df[col].dtype, pd.StringDtype):
                    df[col] = df[col].astype(object)

    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.write(f'{outpath}04.24_BBKNN_scVI.h5ad')

    print(f'scVI-latent BBKNN Done! - {name}')
