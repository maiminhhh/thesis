print('1. UINMF Scanorama Set up')
print('--------------------------------------')

import scanpy as sc
import pandas as pd
import numpy as np
import os
os.environ["OMP_NUM_THREADS"] = os.environ.get("N_CPU", "32")
os.environ["OPENBLAS_NUM_THREADS"] = os.environ.get("N_CPU", "32")
os.environ["MKL_NUM_THREADS"] = os.environ.get("N_CPU", "32")
os.environ["NUMEXPR_NUM_THREADS"] = os.environ.get("N_CPU", "32")
# adata = sc.read_h5ad('data/03_concatenated_data.h5ad')
# H_raw = pd.read_csv('data/04_integrated_data/full_mouse/uinmf_H_raw.csv', index_col = 0)

# print('2. UINMF Scanorama Preprocessing')
# print('--------------------------------------')
# adata.obsm['X_uinmf_H'] = H_raw.loc[adata.obs_names].to_numpy(dtype=float)  

# print('3. Run UINMF Scanorama')
# print('--------------------------------------')
# sc.external.pp.scanorama_integrate(adata, key = 'orig.ident', basis = 'X_uinmf_H', adjusted_basis = 'X_scanorama_uinmf')

# print('4. UINMF Scanorama Postprocessing')
# print('--------------------------------------')
# sc.pp.neighbors(adata, use_rep="X_scanorama_uinmf", n_pcs=None)
# sc.tl.umap(adata, init_pos = 'random')
# adata.obsm['X_umap_scanorama_uinmf'] = adata.obsm['X_umap'].copy()

# print('5. UINMF Scanorama Save data')
# print('--------------------------------------')
# adata.write('data/04_integrated_data/full_mouse/04.15_Scanorama_UINMF.h5ad')

# print('UINMF Scanorama Done!')

# datasets (datasets.py): all by default, or e.g. `python 04.15_integrate_uinmf_scanorama.py genes8` to run one dataset (datasets can run in parallel)
from datasets import selected
genes = selected()
for g in genes:
    adata = sc.read_h5ad(f'data/{g}/03_concatenated_data.h5ad')
    H_raw = pd.read_csv(f'data/{g}/uinmf_H_raw.csv', index_col = 0)

    # same rule as the R UINMF helpers: drop H rows that are non-finite or all zero
    H_raw = H_raw[np.isfinite(H_raw).all(axis=1) & (H_raw.abs().sum(axis=1) > 0)]

    # Cells present in both, kept in adata's order
    common = adata.obs_names.intersection(H_raw.index, sort=False)
    print(f"{g}: keeping {len(common)} / {adata.n_obs} cells "
        f"({adata.n_obs - len(common)} dropped)")

    adata = adata[common].copy()   
    H_raw = H_raw.loc[adata.obs_names]

    adata.obsm['X_uinmf_H'] = H_raw.to_numpy(dtype=float)
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
    adata.write(f'data/{g}/04.15_Scanorama_UINMF.h5ad')

    print('UINMF Scanorama Done!')
