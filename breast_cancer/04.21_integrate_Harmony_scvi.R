print('1. scVI-latent Harmony Set up')
print('--------------------------------------')
# Harmony on the scVI latent (04.18, trained without the species key) instead of PCA
# (../EBCC_PIPELINE.md 6b; adapted from bone_marrow to this tissue's datasets, MERFISH assay and run_umap)
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.21_integrate_Harmony_scvi.R genes8` to run one dataset
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {
data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')
Z = read_scvi_latent(out_path)

print(str_glue('2. Run scVI-latent Harmony - {i}'))
print('--------------------------------------')
data_harmony = integrate_harmony(data, H = Z, W = matrix(numeric(0), 0, 0)) |>
  rename_reduction('uinmf_H', 'scvi_latent')

print(str_glue('3. scVI-latent Harmony Postprocessing - {i}'))
print('--------------------------------------')
data_harmony = data_harmony |>
  run_umap(i, reduction = 'harmony', dims = 1:min(20, ncol(data_harmony[['harmony']])),
          reduction.name = 'umap_harmony_scvi')

print(str_glue('4. scVI-latent Harmony Save data - {i}'))
print('--------------------------------------')
SaveSeuratRds(data_harmony, str_glue('{out_path}04.21_Harmony_scVI.rds'))
as.anndata(data_harmony, file_path = out_path, file_name = '04.21_Harmony_scVI.h5ad',
           assay = 'MERFISH', main_layer = 'counts', other_layers = c('data', 'scale.data'))

print(str_glue('scVI-latent Harmony Done! - {i}'))}
