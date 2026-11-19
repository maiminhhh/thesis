print('1. Scanorama Set up')
print('--------------------------------------')
import os
os.environ.setdefault("OMP_NUM_THREADS", "32")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "32")
os.environ.setdefault("MKL_NUM_THREADS", "32")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "32")
import scanpy as sc
from integration_datasets import DATASETS, tidy_for_write

for name, d in DATASETS.items():
    print(f'===== {name} =====')
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")

    print('2. Scanorama Preprocessing')
    print('--------------------------------------')
    adata.X = adata.layers['data_Xenium'].copy()

    print('3. Run Scanorama')
    print('--------------------------------------')
    sc.pp.pca(adata, layer="scale.data_Xenium", zero_center=False, svd_solver="arpack",
              n_comps=min(50, adata.n_vars - 1))
    sc.external.pp.scanorama_integrate(adata, key = 'orig.ident', basis = 'X_pca', adjusted_basis = 'X_scanorama')

    print('4. Scanorama Postprocessing')
    sc.pp.neighbors(adata, use_rep="X_scanorama", n_pcs=None)
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_scanorama'] = adata.obsm['X_umap'].copy()

    print('5. Scanorama Save data')
    print('--------------------------------------')
    tidy_for_write(adata).write_h5ad(f"{d['out']}04.05_Scanorama.h5ad")

    print('Scanorama Done!')
