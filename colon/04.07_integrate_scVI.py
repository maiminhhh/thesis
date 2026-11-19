print('1. scVI Set up')
print('--------------------------------------')
import scanpy as sc
import scvi
from integration_datasets import DATASETS, tidy_for_write

for name, d in DATASETS.items():
    print(f'===== {name} =====')
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")

    print('2. Run scVI')
    print('--------------------------------------')
    # X holds raw counts (as.anndata main_layer = 'counts')
    adata.layers['counts'] = adata.X
    scvi.model.SCVI.setup_anndata(adata, layer = 'counts', batch_key = 'orig.ident')
    vae = scvi.model.SCVI(adata)
    vae.train()
    adata.obsm['X_scvi'] = vae.get_latent_representation()

    print('4. scVI Postprocessing')
    sc.pp.neighbors(adata, use_rep="X_scvi", n_pcs=None)
    sc.tl.umap(adata, init_pos = 'random')
    adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

    print('3. scVI Save data')
    print('--------------------------------------')
    tidy_for_write(adata).write(f"{d['out']}04.07_scVI.h5ad")

    print('scVI Done!')
