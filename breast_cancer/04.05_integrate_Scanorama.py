# print('1. Scanorama Set up')
# print('--------------------------------------')
# import os
# os.environ["OMP_NUM_THREADS"] = "64"
# os.environ["OPENBLAS_NUM_THREADS"] = "64"
# os.environ["MKL_NUM_THREADS"] = "64"
# os.environ["NUMEXPR_NUM_THREADS"] = "64"
# import click
# import matplotlib.pyplot as plt
# import scanpy as sc
# adata = sc.read_h5ad('data/03_concatenated_data.h5ad')
# adata.X = adata.layers['data_MERFISH'].copy()

# print('2. Scanorama Preprocessing')
# print('--------------------------------------')
# adata.X = adata.layers['data_MERFISH'].copy()

# print('3. Run Scanorama')
# print('--------------------------------------')
# sc.external.pp.scanorama_integrate(adata, key = 'orig.ident', basis = 'X_pca_raw', adjusted_basis = 'X_scanorama')

# print('4. Scanorama Postprocessing')
# sc.pp.neighbors(adata, use_rep="X_scanorama", n_pcs=None)
# sc.tl.umap(adata, init_pos = 'random')
# adata.obsm['X_umap_scanorama'] = adata.obsm['X_umap'].copy()

# print('5. Scanorama Save data')
# print('--------------------------------------')
# adata.write('data/04_integrated_data/full_mouse/04.05_Scanorama.h5ad')

# print('Scanorama Done!')

print('1. Scanorama Set up')
print('--------------------------------------')
import os
import pandas as pd
os.environ["OMP_NUM_THREADS"] = os.environ.get("N_CPU", "32")
os.environ["OPENBLAS_NUM_THREADS"] = os.environ.get("N_CPU", "32")
os.environ["MKL_NUM_THREADS"] = os.environ.get("N_CPU", "32")
os.environ["NUMEXPR_NUM_THREADS"] = os.environ.get("N_CPU", "32")
import scanpy as sc

# datasets (datasets.py): all by default, or e.g. `python 04.05_integrate_Scanorama.py genes8` to run one dataset (datasets can run in parallel)
from datasets import selected
genes_id = selected()

for i in genes_id:
    adata = sc.read_h5ad(f'data/{i}/03_concatenated_data.h5ad')

    print('2. Scanorama Preprocessing')
    print('--------------------------------------')
    adata.X = adata.layers['data_MERFISH'].copy()

    print('3. Run Scanorama')
    print('--------------------------------------')
    # up to 20 PCs, as for the R methods
    sc.pp.pca(adata, layer="scale.data_MERFISH", n_comps=min(20, adata.n_vars - 1), zero_center=False, svd_solver="arpack")
    sc.external.pp.scanorama_integrate(adata, key = 'orig.ident', basis = 'X_pca', adjusted_basis = 'X_scanorama')

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
    adata.obs_names = pd.Index(
        adata.obs_names.astype(str), dtype=object
    )
    adata.var_names = pd.Index(
        adata.var_names.astype(str), dtype=object
    )
    adata.write_h5ad(f'data/{i}/04.05_Scanorama.h5ad')

    print('Scanorama Done!')