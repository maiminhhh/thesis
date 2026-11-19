print('1. scVI-latent RPCA Set up')
print('--------------------------------------')
# Seurat RPCA anchors on the scVI latent (04.18, trained without the species key) instead of PCA
# (../EBCC_PIPELINE.md 6b; adapted from bone_marrow to this tissue's datasets, MERFISH assay and run_umap)
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.23_integrate_RPCA_scvi.R genes8` to run one dataset
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {
data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')
Z = read_scvi_latent(out_path)
# IntegrateEmbeddings carries the input reduction's gene loadings along (it does not
# use them), and needs one column per dimension: the latent has none, so zeros
W = matrix(0, nrow(data), ncol(Z$hsapiens), dimnames = list(rownames(data), NULL))

print(str_glue('2. Run scVI-latent RPCA - {i}'))
print('--------------------------------------')
data_rpca = integrate_rpca(data, H = Z, W = W, half_dims = FALSE) |>
  rename_reduction('rpca_uinmf', 'rpca_scvi') |>
  rename_reduction('uinmf_H', 'scvi_latent')

print(str_glue('3. scVI-latent RPCA Postprocessing - {i}'))
print('--------------------------------------')
data_rpca = data_rpca |>
  run_umap(i, reduction = 'rpca_scvi', dims = 1:min(20, ncol(data_rpca[['rpca_scvi']])),
          reduction.name = 'umap_rpca_scvi')

print(str_glue('4. scVI-latent RPCA Save data - {i}'))
print('--------------------------------------')
SaveSeuratRds(data_rpca, str_glue('{out_path}04.23_RPCA_scVI.rds'))
as.anndata(data_rpca, file_path = out_path, file_name = '04.23_RPCA_scVI.h5ad',
           assay = 'MERFISH', main_layer = 'data', other_layers = 'scale.data')

print(str_glue('scVI-latent RPCA Done! - {i}'))}
