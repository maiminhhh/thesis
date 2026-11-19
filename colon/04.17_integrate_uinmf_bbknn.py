print('1. UINMF BBKNN Setup')
import sys
import sklearn.metrics
import sklearn.neighbors
sklearn.neighbors.DistanceMetric = sklearn.metrics.DistanceMetric
sys.modules['sklearn.neighbors.DistanceMetric'] = sklearn.metrics.DistanceMetric
import scanpy as sc
import bbknn
from integration_datasets import DATASETS, load_uinmf_H, tidy_for_write

for name, d in DATASETS.items():
    print(f'===== {name} =====')
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")
    # H comes from 04.11 and may cover fewer cells: keep only cells with a usable
    # (finite, non-zero) H row, the same rule as the R UINMF scripts
    adata = load_uinmf_H(adata, d['out'])

    print('2. Run UINMF BBKNN')
    bbknn.bbknn(adata, batch_key='orig.ident', use_rep = 'X_uinmf_H')

    print('3. BBKNN Postprocessing')
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_bbknn_uinmf'] = adata.obsm['X_umap']

    print('4. UINMF BBKNN Save data')
    tidy_for_write(adata).write(f"{d['out']}04.17_BBKNN_UINMF.h5ad")

    print('UINMF BBKNN Done!')
