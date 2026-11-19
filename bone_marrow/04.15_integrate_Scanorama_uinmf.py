print('1. UINMF Scanorama Set up')
print('--------------------------------------')
import scanpy as sc
import pandas as pd
import numpy as np
import os
os.environ.setdefault("OMP_NUM_THREADS", "32")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "32")
os.environ.setdefault("MKL_NUM_THREADS", "32")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "32")
from datasets import selected

for name, d in selected():
    outpath = d['out']
    adata = sc.read_h5ad(d['concat'])
    H_raw = pd.read_csv(f'{outpath}uinmf_H_raw.csv', index_col = 0)

    # same cells as the R UINMF methods: in the data, with a finite, non-zero H row
    usable = H_raw.index[np.isfinite(H_raw).all(axis=1) & (H_raw.abs().sum(axis=1) > 0)]
    common = adata.obs_names.intersection(usable, sort=False)
    print(f"{name}: keeping {len(common)} / {adata.n_obs} cells "
        f"({adata.n_obs - len(common)} dropped)")

    adata = adata[common].copy()
    H_raw = H_raw.loc[adata.obs_names]

    adata.obsm['X_uinmf_H'] = H_raw.to_numpy(dtype=float)
    adata.obsm['X_pca_raw'] = adata.obsm['X_pca'].copy()
    adata.obsm['X_umap_raw'] = adata.obsm['X_umap'].copy()

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

    adata.write(f'{outpath}04.15_Scanorama_UINMF.h5ad')

    print(f'UINMF Scanorama Done! - {name}')
