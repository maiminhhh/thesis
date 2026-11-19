# print('1. scVI Set up')
# print('--------------------------------------')
# import matplotlib.pyplot as plt
# import scanpy as sc
# import scvi
# adata = sc.read_h5ad('data/03_concatenated_data.h5ad')

# print('2. Run scVI')
# print('--------------------------------------')
# scvi.model.SCVI.setup_anndata(adata, layer = None, batch_key = 'orig.ident')
# vae = scvi.model.SCVI(adata)
# vae.train()
# vae.history['elbo_train'].plot()
# adata.obsm['X_scvi'] = vae.get_latent_representation()

# print('4. scVI Postprocessing')
# sc.pp.neighbors(adata, use_rep="X_scvi", n_pcs=None)
# sc.tl.umap(adata, init_pos = 'random')
# adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

# print('3. scVI Save data')
# print('--------------------------------------')
# adata.write('data/04_integrated_data/full_mouse/04.07_scVI.h5ad')

# print('scVI Done!')

print('1. scVI Set up')
print('--------------------------------------')
import matplotlib.pyplot as plt
import scanpy as sc
import scvi
import pandas as pd

# datasets (datasets.py): all by default, or e.g. `python 04.07_integrate_scVI.py genes8` to run one dataset (datasets can run in parallel)
from datasets import selected
genes_id = selected()
for i in genes_id:
    adata = sc.read_h5ad(f'data/{i}/03_concatenated_data.h5ad')

    print('2. Run scVI')
    print('--------------------------------------')
    scvi.model.SCVI.setup_anndata(adata, layer = None, batch_key = 'orig.ident')
    vae = scvi.model.SCVI(adata)
    vae.train()
    vae.history['elbo_train'].plot()
    adata.obsm['X_scvi'] = vae.get_latent_representation()
    # 04.08 (scANVI) starts from this model instead of training its own
    vae.save(f'data/{i}/04.07_scVI_model', overwrite=True)

    print('4. scVI Postprocessing')
    sc.pp.neighbors(adata, use_rep="X_scvi", n_pcs=None)
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

    print('3. scVI Save data')
    print('--------------------------------------')
    for df in (adata.obs, adata.var):
        for col in df.columns:
            if isinstance(df[col].dtype, pd.StringDtype):
                df[col] = df[col].astype(object)
    adata.obs_names = pd.Index(adata.obs_names.astype(str), dtype=object)
    adata.var_names = pd.Index(adata.var_names.astype(str), dtype=object)
    adata.write(f'data/{i}/04.07_scVI.h5ad')

    print('scVI Done!')