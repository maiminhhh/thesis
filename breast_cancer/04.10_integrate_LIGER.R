print('1. LIGER Set up')
print('--------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 32)
# data = readRDS('data/03_concatenated_data.rds')
# out_path = 'data/04_integrated_data/full_mouse/'


# print('2. Run LIGER')
# print('--------------------------------------')
# data_liger = integrate_liger(data = data)

# print('3. LIGER Postprocessing')
# data_liger = data_liger |> 
#   RunUMAP(reduction = 'inmfNorm', dims = 1:20,
#           reduction.name = 'umap_liger')

# print('4. LIGER Save data')
# print('--------------------------------------')
# SaveSeuratRds(data_liger, str_glue('{out_path}04.10_LIGER.rds'))
# #write_h5ad(data_liger, str_glue('{out_path}04.10_LIGER.h5ad'))
# as.anndata(data_liger, file_path = out_path, file_name = '04.10_LIGER.h5ad', 
#            assay = 'MERFISH', main_layer = 'counts', other_layers = c('data', 'scale.data'))

# print('LIGER Done!')
# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.10_integrate_LIGER.R genes8` to run one dataset (datasets can run in parallel)
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {
data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')


print('2. Run LIGER')
print('--------------------------------------')
data_liger = integrate_liger(data = data)

print('3. LIGER Postprocessing')
data_liger = data_liger |> 
  run_umap(i, reduction = 'inmfNorm', dims = 1:min(20, ncol(data_liger[['inmfNorm']])),
          reduction.name = 'umap_liger')

print('4. LIGER Save data')
print('--------------------------------------')
SaveSeuratRds(data_liger, str_glue('{out_path}04.10_LIGER.rds'))
#write_h5ad(data_liger, str_glue('{out_path}04.10_LIGER.h5ad'))
as.anndata(data_liger, file_path = out_path, file_name = '04.10_LIGER.h5ad', 
           assay = 'MERFISH', main_layer = 'counts', other_layers = c('data', 'scale.data'))

print('LIGER Done!')}