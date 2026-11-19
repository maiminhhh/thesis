print('1. CCA Set up')
print('--------------------------------------')
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
source('integration_datasets.R')
options(future.globals.maxSize = 20 * 1024^3)

for (d in names(datasets)) {
print(str_glue('===== {d} ====='))
concat = readRDS(str_glue('{datasets[[d]]$concat}.rds'))
out_path = datasets[[d]]$out

print('2. Run CCA')
print('--------------------------------------')
data_cca = integrate_cca(data = concat)
data_cca = data_cca |>
  ScaleData() |>
  RunPCA(reduction.name = 'cca') |>
  RunUMAP(uwot.init = 'pca', reduction = 'cca', dims = 1:min(20, round(dim(data_cca[['cca']])[2]/2)), reduction.name = 'umap_cca')

print('3. CCA Save data')
SaveSeuratRds(data_cca, str_glue('{out_path}04.03_CCA.rds'))
as.anndata(data_cca, file_path = out_path, file_name = '04.03_CCA.h5ad',
           assay = 'integrated', main_layer = 'data', other_layers = 'scale.data')

print('Seurat CCA Done!')}
