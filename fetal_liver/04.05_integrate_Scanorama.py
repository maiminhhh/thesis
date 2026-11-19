print('1. Scanorama Set up')
print('--------------------------------------')
import os
os.environ.setdefault("OMP_NUM_THREADS", "32")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "32")
os.environ.setdefault("MKL_NUM_THREADS", "32")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "32")
import pandas as pd
import scanpy as sc
from integration_datasets import DATASETS
for name, d in DATASETS.items():
    print(f'===== {name} =====')
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")
    outpath = d['out']
    print('2. Scanorama Preprocessing')
    print('--------------------------------------')
    sc.pp.pca(adata, layer = 'scale.data_Xenium')   # scaled log-normalised data, not raw counts
    adata.obsm['X_pca_raw'] = adata.obsm['X_pca']
    print('3. Run Scanorama')
    print('--------------------------------------')
    sc.external.pp.scanorama_integrate(adata, key = 'orig.ident', basis = 'X_pca_raw', adjusted_basis = 'X_scanorama')
    print('4. Scanorama Postprocessing')
    sc.pp.neighbors(adata, use_rep="X_scanorama", n_pcs=None)
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_scanorama'] = adata.obsm['X_umap'].copy()
    print('5. Scanorama Save data')
    print('--------------------------------------')
    for df in (adata.obs, adata.var):
        for col in df.columns:
            if isinstance(df[col].dtype, pd.StringDtype):
                df[col] = df[col].astype(object)
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.write(f'{outpath}04.05_Scanorama.h5ad')
    print('Scanorama Done!')
