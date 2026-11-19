print('1. scVI-latent Harmony setup')
print('--------------------------------------')
# Harmony on the scVI latent (04.18, no species key) instead of PCA
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
source('integration_datasets.R')

for (s in names(datasets)) {
data = readRDS(str_glue('{datasets[[s]]$concat}.rds'))
out_path = datasets[[s]]$out
Z = read_scvi_latent(out_path)

print(str_glue('2. Run scVI-latent Harmony - {s}'))
print('--------------------------------------')
data_harmony = integrate_harmony(data, H = Z, W = matrix(numeric(0), 0, 0)) |>
  rename_reduction('uinmf_H', 'scvi_latent')

print('3. scVI-latent Harmony Postprocessing')
print('--------------------------------------')
data_harmony = data_harmony |>
  RunUMAP(reduction = 'harmony',
          dims = 1:min(20, ncol(data_harmony[['harmony']])),
          reduction.name = 'umap_harmony_scvi')

print('4. scVI-latent Harmony Save data')
print('--------------------------------------')
SaveSeuratRds(data_harmony, str_glue('{out_path}04.21_Harmony_scVI.rds'))
as.anndata(data_harmony, file_path = out_path,
           file_name = '04.21_Harmony_scVI.h5ad',
           assay = 'Xenium', main_layer = 'counts', other_layers = c('data', 'scale.data'))

print(str_glue('scVI-latent Harmony Done! - {s}'))
}
