print('1. UINMF BBKNN Setup')
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
    H_raw = pd.read_csv(f'{outpath}uinmf_H_raw.csv', index_col = 0)

    # same cells as the R UINMF methods: in the data, with a finite, non-zero H row
    usable = H_raw.index[np.isfinite(H_raw).all(axis=1) & (H_raw.abs().sum(axis=1) > 0)]
    common = adata.obs_names.intersection(usable, sort=False)
    print(f"{name}: keeping {len(common)} / {adata.n_obs} cells "
        f"({adata.n_obs - len(common)} dropped)")
    adata = adata[common].copy()
    adata.obsm['X_uinmf_H'] = H_raw.loc[adata.obs_names].to_numpy(dtype=float)

    print('Preprocessing')
    adata.obsm['X_raw_pca'] = adata.obsm['X_pca'].copy()
    adata.obsm['X_raw_umap'] = adata.obsm['X_umap'].copy()

    print('2. Run UINMF BBKNN')
    bbknn.bbknn(adata, batch_key='orig.ident', use_rep = 'X_uinmf_H')

    print('3. BBKNN Postprocessing')
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_bbknn_uinmf'] = adata.obsm['X_umap']

    print('4. UINMF BBKNN Save data')
    for df in (adata.obs, adata.var):
            for col in df.columns:
                if isinstance(df[col].dtype, pd.StringDtype):
                    df[col] = df[col].astype(object)

    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.write(f'{outpath}04.17_BBKNN_UINMF.h5ad')

    print(f'UINMF BBKNN Done! - {name}')
