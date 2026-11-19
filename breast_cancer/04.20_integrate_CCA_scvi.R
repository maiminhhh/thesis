print('1. scVI-latent CCA Set up')
print('--------------------------------------')
# Seurat CCA anchors on the scVI latent (04.18, trained without the species key) instead of PCA
# (../EBCC_PIPELINE.md 6b; adapted from bone_marrow to this tissue's datasets, MERFISH assay and run_umap)
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.20_integrate_CCA_scvi.R genes8` to run one dataset
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {
data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')
Z = read_scvi_latent(out_path)
# IntegrateEmbeddings carries the input reduction's gene loadings along (it does not
# use them), and needs one column per dimension: the latent has none, so zeros
W = matrix(0, nrow(data), ncol(Z$hsapiens), dimnames = list(rownames(data), NULL))

print(str_glue('2. Run scVI-latent CCA - {i}'))
print('--------------------------------------')
data_cca = integrate_cca(data = data, H = Z, W = W, half_dims = FALSE) |>
  rename_reduction('cca_uinmf', 'cca_scvi') |>
  rename_reduction('uinmf_H', 'scvi_latent')
# replace non-finite values in the integrated embedding before UMAP
emb = Embeddings(data_cca, 'cca_scvi')
emb[!is.finite(emb)] = 0
data_cca[['cca_scvi']] = CreateDimReducObject(embeddings = emb, key = 'ccascvi_', assay = DefaultAssay(data_cca))

print(str_glue('3. scVI-latent CCA Postprocessing - {i}'))
print('--------------------------------------')
data_cca = data_cca |>
  run_umap(i, reduction = 'cca_scvi', dims = 1:min(20, ncol(data_cca[['cca_scvi']])),
          reduction.name = 'umap_cca_scvi')

print(str_glue('4. scVI-latent CCA Save data - {i}'))
print('--------------------------------------')
SaveSeuratRds(data_cca, str_glue('{out_path}04.20_CCA_scVI.rds'))
as.anndata(data_cca, file_path = out_path, file_name = '04.20_CCA_scVI.h5ad',
           assay = 'MERFISH', main_layer = 'data', other_layers = 'scale.data')

print(str_glue('scVI-latent CCA Done! - {i}'))}
