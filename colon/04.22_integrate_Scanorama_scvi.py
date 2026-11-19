print('1. scVI-latent Scanorama Set up')
# Scanorama on the scVI latent (04.18, no species key) instead of PCA
print('--------------------------------------')
import scanpy as sc
import pandas as pd
import numpy as np
import os
os.environ.setdefault("OMP_NUM_THREADS", "32")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "32")
os.environ.setdefault("MKL_NUM_THREADS", "32")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "32")
from integration_datasets import DATASETS

for name, d in DATASETS.items():
    outpath = d['out']
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")
    Z = pd.read_csv(f'{outpath}scvi_latent.csv', index_col = 0)
    adata.obsm['X_scvi_latent'] = Z.loc[adata.obs_names].to_numpy(dtype=float)
    adata.obsm['X_pca_raw'] = adata.obsm['X_pca'].copy()
    adata.obsm['X_umap_raw'] = adata.obsm['X_umap'].copy()

    print('3. Run scVI-latent Scanorama')
    print('--------------------------------------')
    sc.external.pp.scanorama_integrate(adata, key = 'orig.ident', basis = 'X_scvi_latent', adjusted_basis = 'X_scanorama_scvi')

    print('4. scVI-latent Scanorama Postprocessing')
    print('--------------------------------------')
    sc.pp.neighbors(adata, use_rep="X_scanorama_scvi", n_pcs=None)
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_scanorama_scvi'] = adata.obsm['X_umap'].copy()

    print('5. scVI-latent Scanorama Save data')
    print('--------------------------------------')
    for df in (adata.obs, adata.var):
            for col in df.columns:
                if isinstance(df[col].dtype, pd.StringDtype):
                    df[col] = df[col].astype(object)

    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)

    adata.write(f'{outpath}04.22_Scanorama_scVI.h5ad')

    print(f'scVI-latent Scanorama Done! - {name}')
