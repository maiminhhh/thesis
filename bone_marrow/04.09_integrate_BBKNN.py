print('1. BBKNN Setup')
import sys
import sklearn.metrics
import sklearn.neighbors
sklearn.neighbors.DistanceMetric = sklearn.metrics.DistanceMetric
sys.modules['sklearn.neighbors.DistanceMetric'] = sklearn.metrics.DistanceMetric
import scanpy as sc
import bbknn
import pandas as pd
from datasets import selected

for name, d in selected():
    print(f'Load data - {name}:')
    adata = sc.read_h5ad(d['concat'])

    print('Preprocessing')
    sc.pp.pca(adata, layer = 'scale.data_Xenium')   # scaled log-normalised data, not raw counts

    print('2. Run BBKNN')
    bbknn.bbknn(adata, batch_key='orig.ident')

    print('3. BBKNN Postprocessing')
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_bbknn'] = adata.obsm['X_umap']

    print('4. BBKNN Save data')
    for df in (adata.obs, adata.var):
        for col in df.columns:
            if isinstance(df[col].dtype, pd.StringDtype):
                df[col] = df[col].astype(object)

    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.write(f"{d['out']}04.09_BBKNN.h5ad")

    print(f'BBKNN Done! - {name}')
