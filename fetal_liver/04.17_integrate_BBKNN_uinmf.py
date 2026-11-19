print('1. UINMF BBKNN Setup')
import sys
import sklearn.metrics
import sklearn.neighbors
sklearn.neighbors.DistanceMetric = sklearn.metrics.DistanceMetric
sys.modules['sklearn.neighbors.DistanceMetric'] = sklearn.metrics.DistanceMetric
import scanpy as sc
import bbknn
import pandas as pd
from integration_datasets import DATASETS, load_uinmf_H
for name, d in DATASETS.items():
    print(f'===== {name} =====')
    print('Load data:')
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")
    print('Preprocessing')
    # keeps only cells with a usable (finite, non-zero) H row, as the R UINMF scripts do
    adata = load_uinmf_H(adata, d['out'])
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
    adata.write(f"{d['out']}04.17_BBKNN_UINMF.h5ad")
    print('UINMF BBKNN Done!')
