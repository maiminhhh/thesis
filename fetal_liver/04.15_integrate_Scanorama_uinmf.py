print('1. UINMF Scanorama Set up')
print('--------------------------------------')
import scanpy as sc
import pandas as pd
import os
os.environ.setdefault("OMP_NUM_THREADS", "32")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "32")
os.environ.setdefault("MKL_NUM_THREADS", "32")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "32")
from integration_datasets import DATASETS, load_uinmf_H
for name, d in DATASETS.items():
    print(f'===== {name} =====')
    print('2. UINMF Scanorama Preprocessing')
    print('--------------------------------------')
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")
    # keeps only cells with a usable (finite, non-zero) H row, as the R UINMF scripts do
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
    for df in (adata.obs, adata.var):
        for col in df.columns:
            if isinstance(df[col].dtype, pd.StringDtype):
                df[col] = df[col].astype(object)
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.write(f"{d['out']}04.15_Scanorama_UINMF.h5ad")
    print('UINMF Scanorama Done!')
