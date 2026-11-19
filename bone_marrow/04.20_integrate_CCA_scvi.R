print('1. scVI-latent CCA Set up')
print('--------------------------------------')
# Seurat CCA anchors on the scVI latent (04.18, no species key) instead of PCA
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
source('integration_sets.R')

for (s in run_sets) {
concat = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
Z = read_scvi_latent(out_path)
# IntegrateEmbeddings carries the input reduction's gene loadings along (it does not
# use them), and needs one column per dimension: the latent has none, so zeros
W = matrix(0, nrow(concat), ncol(Z$hsapiens), dimnames = list(rownames(concat), NULL))

print(str_glue('2. Run scVI-latent CCA - {s}'))
print('--------------------------------------')
data_cca = integrate_cca(data = concat, H = Z, W = W) |>
  rename_reduction('cca_uinmf', 'cca_scvi') |>
  rename_reduction('uinmf_H', 'scvi_latent')

print('3. scVI-latent CCA Postprocessing')
print('--------------------------------------')
# replace non-finite values in the integrated embedding before UMAP
emb = Embeddings(data_cca, 'cca_scvi')
emb[!is.finite(emb)] = 0
data_cca[['cca_scvi']] = CreateDimReducObject(embeddings = emb, key = 'ccascvi_',
                                              assay = DefaultAssay(data_cca))

data_cca = data_cca |>
  RunUMAP(reduction = 'cca_scvi', dims = 1:min(20, round(dim(data_cca[['cca_scvi']])[2]/2)),
          reduction.name = 'umap_cca_scvi')

print('4. scVI-latent CCA Save data')
SaveSeuratRds(data_cca, str_glue('{out_path}04.20_CCA_scVI.rds'))
as.anndata(data_cca, file_path = out_path, file_name = '04.20_CCA_scVI.h5ad',
           assay = 'Xenium', main_layer = 'data', other_layers = 'scale.data')

print(str_glue('scVI-latent Seurat CCA Done! - {s}'))
}
