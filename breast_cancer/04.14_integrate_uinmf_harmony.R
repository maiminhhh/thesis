
print('1. UINMF Harmony setup')
print('--------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 32)
# data = readRDS('data/03_concatenated_data.rds')
# out_path = 'data/04_integrated_data/full_mouse/'

# H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
# W = readRDS(str_glue('{out_path}uinmf_W.rds'))

# print('2. Run UINMF Harmony')
# print('--------------------------------------')
# data_harmony = integrate_harmony(data, H = H, W = W)

# print('3. UINMF Harmony Postprocessing')
# print('--------------------------------------')
# data_harmony = data_harmony |> 
#   RunUMAP(reduction = 'harmony', 
#           dims = 1:20,
#           reduction.name = 'umap_harmony_uinmf')

# print('3. UINMF Harmony Save data')
# print('--------------------------------------')
# SaveSeuratRds(data_harmony, str_glue('{out_path}04.14_Harmony_UINMF.rds'))
# #write_h5ad(data_harmony, str_glue('{out_path}04.14_Harmony_UINMF.h5ad'))
# as.anndata(data_harmony, file_path = out_path, 
#            file_name = '04.14_Harmony_UINMF.h5ad',
#            assay = 'MERFISH', main_layer = 'counts', other_layers = c('data', 'scale.data'))

# print('UINMF Harmony Done!')

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.14_integrate_uinmf_harmony.R genes8` to run one dataset (datasets can run in parallel)
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {
data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')

H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print('2. Run UINMF Harmony')
print('--------------------------------------')
data_harmony = integrate_harmony(data, H = H, W = W)

print('3. UINMF Harmony Postprocessing')
print('--------------------------------------')
data_harmony = data_harmony |> 
  run_umap(i, reduction = 'harmony', 
          dims = 1:min(20, ncol(data_harmony[['harmony']])),
          reduction.name = 'umap_harmony_uinmf')

print('3. UINMF Harmony Save data')
print('--------------------------------------')
SaveSeuratRds(data_harmony, str_glue('{out_path}04.14_Harmony_UINMF.rds'))
#write_h5ad(data_harmony, str_glue('{out_path}04.14_Harmony_UINMF.h5ad'))
as.anndata(data_harmony, file_path = out_path, 
           file_name = '04.14_Harmony_UINMF.h5ad',
           assay = 'MERFISH', main_layer = 'counts', other_layers = c('data', 'scale.data'))

print('UINMF Harmony Done!')}