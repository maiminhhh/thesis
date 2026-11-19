print('1. scVI latent (no species key) Set up')
print('--------------------------------------')
# scVI trained WITHOUT a batch key: a denoised, non-linear counterpart of PCA in which the
# species are still separate. 04.19-04.24 run FastMNN, CCA, Harmony, Scanorama, RPCA and
# BBKNN on this latent instead of PCA (as 04.12-04.17 do on the UINMF H matrix).
# 20 latent dimensions to match the 20 UINMF factors.
import scanpy as sc
import scvi
import pandas as pd
from integration_datasets import DATASETS

N_LATENT = 20
scvi.settings.seed = 1

for name, d in DATASETS.items():
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")

    print(f'2. Run scVI (no species key) - {name}')
    print('--------------------------------------')
    adata.layers['counts'] = adata.X
    scvi.model.SCVI.setup_anndata(adata, layer = 'counts')
    vae = scvi.model.SCVI(adata, n_latent = N_LATENT)
    vae.train()
    adata.obsm['X_scvi_latent'] = vae.get_latent_representation()
    # cells x latent dimensions, read by the R and Python methods built on it
    pd.DataFrame(adata.obsm['X_scvi_latent'], index = adata.obs_names,
                 columns = [f'scvi_{i + 1}' for i in range(N_LATENT)]).to_csv(f"{d['out']}scvi_latent.csv")

    print('3. scVI latent Postprocessing')
    sc.pp.neighbors(adata, use_rep = 'X_scvi_latent', n_pcs = None)
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_scvi_latent'] = adata.obsm['X_umap'].copy()

    print('4. scVI latent Save data')
    print('--------------------------------------')
    for df in (adata.obs, adata.var):
        for col in df.columns:
            if isinstance(df[col].dtype, pd.StringDtype):
                df[col] = df[col].astype(object)

    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.write(f"{d['out']}04.18_scVI_latent.h5ad")

    print(f'scVI latent Done! - {name}')
