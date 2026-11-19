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

# print('Load data:')
# adata = sc.read_h5ad('data/03_concatenated_data.h5ad')
# H_raw = pd.read_csv('data/04_integrated_data/full_mouse/uinmf_H_raw.csv', index_col = 0)
# adata.obsm['X_uinmf_H'] = H_raw.loc[adata.obs_names].to_numpy(dtype=float)  

# print('Preprocessing')

# print('2. Run UINMF BBKNN')
# bbknn.bbknn(adata, batch_key='orig.ident', use_rep = 'X_uinmf_H')

# print('3. BBKNN Postprocessing')
# sc.tl.umap(adata, init_pos = 'random')
# adata.obsm['X_umap_bbknn_uinmf'] = adata.obsm['X_umap']

# print('4. UINMF BBKNN Save data')
# adata.write('data/04_integrated_data/full_mouse/04.17_BBKNN_UINMF.h5ad')

# print('UINMF BBKNN Done!')

# datasets (datasets.py): all by default, or e.g. `python 04.17_integrate_uinmf_bbknn.py genes8` to run one dataset (datasets can run in parallel)
from datasets import selected
genes = selected()
for g in genes:
    print('Load data:')
    adata = sc.read_h5ad(f'data/{g}/03_concatenated_data.h5ad')
    H_raw = pd.read_csv(f'data/{g}/uinmf_H_raw.csv', index_col = 0)
    # Reorder H to match AnnData
    # keep only cells with an H row (reindex would give NaN rows to the rest);
    # same rule as the R UINMF helpers: drop H rows that are non-finite or all zero
    H_raw = H_raw[np.isfinite(H_raw).all(axis=1) & (H_raw.abs().sum(axis=1) > 0)]
    common = adata.obs_names.intersection(H_raw.index, sort=False)
    print(f"{g}: keeping {len(common)} / {adata.n_obs} cells")
    adata = adata[common].copy()
    print('2. UINMF Scanorama Preprocessing')
    print('--------------------------------------')
    adata.obsm['X_uinmf_H'] = H_raw.loc[adata.obs_names].to_numpy(dtype=float)  

    print('Preprocessing')

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
    adata.write(f'data/{g}/04.17_BBKNN_UINMF.h5ad')

    print('UINMF BBKNN Done!')