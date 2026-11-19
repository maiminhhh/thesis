print('1. UINMF Scanorama Set up')
print('--------------------------------------')
import os
os.environ.setdefault("OMP_NUM_THREADS", "32")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "32")
os.environ.setdefault("MKL_NUM_THREADS", "32")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "32")
import scanpy as sc
from integration_datasets import DATASETS, load_uinmf_H, tidy_for_write

for name, d in DATASETS.items():
    print(f'===== {name} =====')
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")
    # H comes from 04.11 and may cover fewer cells: keep only cells with a usable
    # (finite, non-zero) H row, the same rule as the R UINMF scripts
    adata = load_uinmf_H(adata, d['out'])

    print('3. Run UINMF Scanorama')
    print('--------------------------------------')
    sc.external.pp.scanorama_integrate(adata, key = 'orig.ident', basis = 'X_uinmf_H', adjusted_basis = 'X_scanorama_uinmf')

    print('4. UINMF Scanorama Postprocessing')
    print('--------------------------------------')
    sc.pp.neighbors(adata, use_rep="X_scanorama_uinmf", n_pcs=None)
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_scanorama_uinmf'] = adata.obsm['X_umap'].copy()

    print('5. UINMF Scanorama Save data')
    print('--------------------------------------')
    tidy_for_write(adata).write(f"{d['out']}04.15_Scanorama_UINMF.h5ad")

    print('UINMF Scanorama Done!')
