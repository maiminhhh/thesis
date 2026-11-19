print('1. BBKNN Setup')
import sys
import sklearn.metrics
import sklearn.neighbors
sklearn.neighbors.DistanceMetric = sklearn.metrics.DistanceMetric
sys.modules['sklearn.neighbors.DistanceMetric'] = sklearn.metrics.DistanceMetric
import scanpy as sc
import bbknn
import pandas as pd

# print('Load data:')
# adata = sc.read_h5ad('data/03_concatenated_data.h5ad')

# print('Preprocessing')
# adata.obsm['X_pca'] = adata.obsm['X_pca_raw'].copy()

# print('2. Run BBKNN')
# bbknn.bbknn(adata, batch_key='orig.ident')

# print('3. BBKNN Postprocessing')
# sc.tl.umap(adata, init_pos = 'random')
# adata.obsm['X_umap_bbknn'] = adata.obsm['X_umap']

# print('4. BBKNN Save data')
# adata.write('data/04_integrated_data/full_mouse/04.09_BBKNN.h5ad')

# print('BBKNN Done!')

# datasets (datasets.py): all by default, or e.g. `python 04.09_integrate_BBKNN.py genes8` to run one dataset (datasets can run in parallel)
from datasets import selected
genes_id = selected()
for i in genes_id:
    print('Load data:')
    adata = sc.read_h5ad(f'data/{i}/03_concatenated_data.h5ad')

    print('Preprocessing')
    # up to 20 PCs, as for the R methods
    sc.pp.pca(adata, layer="scale.data_MERFISH", n_comps=min(20, adata.n_vars - 1), zero_center=False, svd_solver="arpack")

    print('2. Run BBKNN')
    bbknn.bbknn(adata, batch_key='orig.ident', n_pcs=adata.obsm['X_pca'].shape[1])

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

    adata.write(f'data/{i}/04.09_BBKNN.h5ad')

    print('BBKNN Done!')