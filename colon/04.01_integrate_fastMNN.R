print('1. FastMNN Set up')
print('--------------------------------------')
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
source('integration_datasets.R')

for (d in names(datasets)) {
print(str_glue('===== {d} ====='))
data = readRDS(str_glue('{datasets[[d]]$concat}.rds'))
out_path = datasets[[d]]$out
print('2. Run FastMNN')
print('--------------------------------------')
data_fastmnn = integrate_fastmnn(data)
data_fastmnn = JoinLayers(data_fastmnn)

data_fastmnn = data_fastmnn |>
  RunUMAP(uwot.init = 'pca', reduction = 'mnn', dims = 1:min(20, round(dim(data_fastmnn[['mnn']])[2]/2)),
          reduction.name = 'umap_mnn')

print('3. FastMNN Save data')
print('--------------------------------------')
SaveSeuratRds(data_fastmnn, str_glue('{out_path}04.01_FastMNN.rds'))

as.anndata(data_fastmnn, file_path = out_path, file_name = '04.01_FastMNN.h5ad',
           assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))

print('FastMNN Done!')}
