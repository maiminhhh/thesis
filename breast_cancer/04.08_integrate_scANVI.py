print('1. scANVI Setup')
print('--------------------------------------')
import os
import scanpy as sc
import scvi
import pandas as pd
# adata = sc.read_h5ad('data/03_concatenated_data.h5ad')

# print('3. Run scANVI')
# print('--------------------------------------')
# scvi.model.SCVI.setup_anndata(adata, batch_key='orig.ident', layer=None)
# scvi_model = scvi.model.SCVI(adata)
# scvi_model.train()

# scanvi_model = scvi.model.SCANVI.from_scvi_model(
#     scvi_model, labels_key='ebcc_label', unlabeled_category='Unknown')
# scanvi_model.train()

# print('4. scANVI Postprocessing')
# print('--------------------------------------')
# adata.obsm['X_scanvi'] = scanvi_model.get_latent_representation()
# adata.obs['scanvi_predicted'] = scanvi_model.predict()
# pred_probs = scanvi_model.predict(soft=True)
# adata.obs['scanvi_confidence'] = pred_probs.max(axis=1).values

# sc.pp.neighbors(adata, use_rep="X_scanvi", n_pcs=None)
# sc.tl.umap(adata, init_pos = 'random')
# adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()

# print('4. scANVI Save data')
# print('--------------------------------------')
# adata.write('data/04_integrated_data/full_mouse/04.08_scANVI.h5ad')

# print('scANVI Done!')

# datasets (datasets.py): all by default, or e.g. `python 04.08_integrate_scANVI.py genes8` to run one dataset (datasets can run in parallel)
from datasets import selected
genes_id = selected()
for i in genes_id:
    adata = sc.read_h5ad(f'data/{i}/03_concatenated_data.h5ad')

    print('3. Run scANVI')
    print('--------------------------------------')
    scvi.model.SCVI.setup_anndata(adata, batch_key='orig.ident', layer=None)
    # reuse the scVI model trained by 04.07 (run it first); train one only if it is missing
    model_dir = f'data/{i}/04.07_scVI_model'
    if os.path.isdir(model_dir):
        print(f'Loading scVI model from {model_dir}')
        scvi_model = scvi.model.SCVI.load(model_dir, adata=adata)
    else:
        print(f'{model_dir} not found: training scVI from scratch')
        scvi_model = scvi.model.SCVI(adata)
        scvi_model.train()

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, labels_key='ebcc_label', unlabeled_category='Unknown')
    scanvi_model.train()

    print('4. scANVI Postprocessing')
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
    adata.obs_names = pd.Index(adata.obs_names.astype(str), dtype=object)
    adata.var_names = pd.Index(adata.var_names.astype(str), dtype=object)
    adata.write(f'data/{i}/04.08_scANVI.h5ad')

    print('scANVI Done!')