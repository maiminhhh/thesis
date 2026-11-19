print('1. RPCA Set up')
print('--------------------------------------')
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'fetal_liver'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
source('integration_datasets.R')

for (d in names(datasets)) {
print(str_glue('===== {d} ====='))
print('2. Run RPCA')
data = readRDS(str_glue('{datasets[[d]]$concat}.rds'))
out_path = datasets[[d]]$out
print('--------------------------------------')
# integrate_rpca already scales the integrated assay and builds the 'rpca' PCA
data_rpca = integrate_rpca(data)

print('3. RPCA Postprocessing')
print('--------------------------------------')
data_rpca = data_rpca |>
  RunUMAP(uwot.init = 'pca', reduction = 'rpca', dims = 1:min(20, round(dim(data_rpca[['rpca']])[2]/2)),
          reduction.name = 'umap_rpca')

print('4. RPCA Save data')
print('--------------------------------------')
SaveSeuratRds(data_rpca, str_glue('{out_path}04.06_RPCA.rds'))
as.anndata(data_rpca, file_path = out_path, file_name = '04.06_RPCA.h5ad',
           assay = 'integrated', main_layer = 'data', other_layers = 'scale.data')

print('Seurat RPCA Done!')}
