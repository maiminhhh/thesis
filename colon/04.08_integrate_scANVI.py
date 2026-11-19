print('1. scANVI Setup')
print('--------------------------------------')
import scanpy as sc
import scvi
from integration_datasets import DATASETS, tidy_for_write

for name, d in DATASETS.items():
    print(f'===== {name} =====')
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")

    print('3. Run scANVI')
    print('--------------------------------------')
    scvi.model.SCVI.setup_anndata(adata, batch_key='orig.ident', layer=None)
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
    tidy_for_write(adata).write(f"{d['out']}04.08_scANVI.h5ad")

    print('scANVI Done!')
