print('1. scVI Set up')
print('--------------------------------------')
import scanpy as sc
import scvi
import pandas as pd
from integration_datasets import DATASETS
for name, d in DATASETS.items():
    print(f'===== {name} =====')
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")
    print('2. Run scVI')
    print('--------------------------------------')
    # X holds raw counts (as.anndata main_layer = 'counts')
    scvi.model.SCVI.setup_anndata(adata, layer = None, batch_key = 'orig.ident')
    vae = scvi.model.SCVI(adata)
    vae.train()
    # saved so 04.08 (scANVI) starts from this model instead of training its own
    vae.save(f"{d['out']}04.07_scVI_model", overwrite = True)
    adata.obsm['X_scvi'] = vae.get_latent_representation()
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
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.write(f"{d['out']}04.07_scVI.h5ad")
    print('scVI Done!')
