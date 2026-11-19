# print('1. Harmony setup')
# print('--------------------------------------')
# source('../library.R')
# source('../source.R')
# param = MulticoreParam(workers = 64)
# data = readRDS('data/03_concatenated_data.rds')
# out_path = 'data/04_integrated_data/full_mouse/'

# print('2. Run Harmony')
# print('--------------------------------------')
# data_harmony = integrate_harmony(data)
# data_harmony = data_harmony |> 
#   RunUMAP(reduction = 'harmony', dims = 1:20,
#           reduction.name = 'umap_harmony')

# print('3. Save data')
# print('--------------------------------------')
# SaveSeuratRds(data_harmony, str_glue('{out_path}04.04_Harmony.rds'))
# #write_h5ad(data_harmony, str_glue('{out_path}04.04_Harmony.h5ad'))
# as.anndata(data_harmony, file_path = out_path, file_name = '04.04_Harmony.h5ad',
#            assay = 'MERFISH', main_layer = 'counts', other_layers = c('data', 'scale.data'))

# print('Harmony Done!')

print('1. Harmony setup')
print('--------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 32)

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.04_integrate_Harmony.R genes8` to run one dataset (datasets can run in parallel)
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {

data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')

print('2. Run Harmony')
print('--------------------------------------')
data_harmony = integrate_harmony(data, max_pcs = 20)
data_harmony = data_harmony |> 
  run_umap(i, reduction = 'harmony', dims = 1:min(20, ncol(data_harmony[['harmony']])),
          reduction.name = 'umap_harmony')

print('3. Save data')
print('--------------------------------------')
SaveSeuratRds(data_harmony, str_glue('{out_path}04.04_Harmony.rds'))
#write_h5ad(data_harmony, str_glue('{out_path}04.04_Harmony.h5ad'))
as.anndata(data_harmony, file_path = out_path, file_name = '04.04_Harmony.h5ad',
           assay = 'MERFISH', main_layer = 'counts', other_layers = c('data', 'scale.data'))

print('Harmony Done!')
}