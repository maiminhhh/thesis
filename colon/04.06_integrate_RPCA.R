print('1. RPCA Set up')
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
# 03 scales both species jointly; RPCA (scale = FALSE in integrate_rpca) expects each species scaled on its own
data = data |>
  ScaleData(features = rownames(data), split.by = 'orig.ident')

print('2. Run RPCA')
print('--------------------------------------')
data_rpca = integrate_rpca(data)

print('3. RPCA Postprocessing')
print('--------------------------------------')
DefaultAssay(data_rpca) = 'integrated'
data_rpca = data_rpca |>
  ScaleData() |>
  RunPCA(reduction.name = 'rpca') |>
  RunUMAP(uwot.init = 'pca', reduction = 'rpca', dims = 1:min(20, round(dim(data_rpca[['rpca']])[2]/2)),
          reduction.name = 'umap_rpca')

print('4. RPCA Save data')
print('--------------------------------------')
SaveSeuratRds(data_rpca, str_glue('{out_path}04.06_RPCA.rds'))
as.anndata(data_rpca, file_path = out_path, file_name = '04.06_RPCA.h5ad',
           assay = 'integrated', main_layer = 'data', other_layers = 'scale.data')

print('Seurat RPCA Done!')}
