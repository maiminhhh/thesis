print('1. LIGER Set up')
print('--------------------------------------')
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'fetal_liver'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
source('integration_datasets.R')

for (d in names(datasets)) {
print(str_glue('===== {d} ====='))
print('2. Run LIGER')
data = readRDS(str_glue('{datasets[[d]]$concat}.rds'))
out_path = datasets[[d]]$out
data_liger = integrate_liger(data = data, n_cores = n_threads)

print('3. LIGER Postprocessing')
data_liger = data_liger |>
  RunUMAP(uwot.init = 'pca', reduction = 'inmfNorm', dims = 1:min(20, round(dim(data_liger[['inmfNorm']])[2]/2)),
          reduction.name = 'umap_liger')

print('4. LIGER Save data')
print('--------------------------------------')
SaveSeuratRds(data_liger, str_glue('{out_path}04.10_LIGER.rds'))
as.anndata(data_liger, file_path = out_path, file_name = '04.10_LIGER.h5ad',
           assay = 'Xenium', main_layer = 'counts', other_layers = c('data', 'scale.data'))

print('LIGER Done!')
}
