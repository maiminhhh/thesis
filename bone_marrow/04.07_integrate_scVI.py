print('1. scVI Set up')
print('--------------------------------------')
import matplotlib.pyplot as plt
import scanpy as sc
import scvi
import pandas as pd
from datasets import selected

for name, d in selected():
    adata = sc.read_h5ad(d['concat'])

    print(f'2. Run scVI - {name}')
    print('--------------------------------------')
    adata.layers['counts'] = adata.X
    scvi.model.SCVI.setup_anndata(adata, layer = 'counts', batch_key = 'orig.ident')
    vae = scvi.model.SCVI(adata)
    vae.train()
    vae.history['elbo_train'].plot()
    adata.obsm['X_scvi'] = vae.get_latent_representation()
    # 04.08 builds scANVI on this model instead of training another scVI
    vae.save(f"{d['out']}04.07_scVI_model", overwrite = True)

    print('3. scVI Postprocessing')
    sc.pp.neighbors(adata, use_rep="X_scvi", n_pcs=None)
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

    print('4. scVI Save data')
    print('--------------------------------------')
    for df in (adata.obs, adata.var):
        for col in df.columns:
            if isinstance(df[col].dtype, pd.StringDtype):
                df[col] = df[col].astype(object)

    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.write(f"{d['out']}04.07_scVI.h5ad")

    print(f'scVI Done! - {name}')
