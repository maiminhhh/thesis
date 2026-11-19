print('1. scANVI Setup')
print('--------------------------------------')
import matplotlib.pyplot as plt
import scanpy as sc
import scvi
import pandas as pd
from datasets import selected

for name, d in selected():
    adata = sc.read_h5ad(d['concat'])

    print(f'2. Run scANVI - {name}')
    print('--------------------------------------')
    adata.layers['counts'] = adata.X
    # the scVI model trained (and evaluated) in 04.07, so it is not trained twice
    scvi_model = scvi.model.SCVI.load(f"{d['out']}04.07_scVI_model", adata = adata)

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, labels_key='ebcc_label', unlabeled_category='Unknown')
    scanvi_model.train()

    print('3. scANVI Postprocessing')
    print('--------------------------------------')
    adata.obsm['X_scanvi'] = scanvi_model.get_latent_representation()
    adata.obs['scanvi_predicted'] = scanvi_model.predict()
    pred_probs = scanvi_model.predict(soft=True)
    adata.obs['scanvi_confidence'] = pred_probs.max(axis=1).values

    sc.pp.neighbors(adata, use_rep="X_scanvi", n_pcs=None)
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()

    print('4. scANVI Save data')
    print('--------------------------------------')
    for df in (adata.obs, adata.var):
        for col in df.columns:
            if isinstance(df[col].dtype, pd.StringDtype):
                df[col] = df[col].astype(object)

    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.write(f"{d['out']}04.08_scANVI.h5ad")

    print(f'scANVI Done! - {name}')
