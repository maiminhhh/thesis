print('1. BBKNN Setup')
import sys
import sklearn.metrics
import sklearn.neighbors
sklearn.neighbors.DistanceMetric = sklearn.metrics.DistanceMetric
sys.modules['sklearn.neighbors.DistanceMetric'] = sklearn.metrics.DistanceMetric
import scanpy as sc
import bbknn
from integration_datasets import DATASETS, tidy_for_write

for name, d in DATASETS.items():
    print(f'===== {name} =====')
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")

    print('Preprocessing')
    # scaled log-normalised data, not raw counts
    sc.pp.pca(adata, layer = 'scale.data_Xenium', n_comps = min(50, adata.n_vars - 1))

    print('2. Run BBKNN')
    bbknn.bbknn(adata, batch_key='orig.ident')

    print('3. BBKNN Postprocessing')
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_bbknn'] = adata.obsm['X_umap']

    print('4. BBKNN Save data')
    tidy_for_write(adata).write(f"{d['out']}04.09_BBKNN.h5ad")

    print('BBKNN Done!')
