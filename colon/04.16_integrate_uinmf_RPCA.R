print('1. UINMF RPCA Set up')
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

# H comes from 04.11 and may cover fewer cells; integrate_rpca keeps only cells
# with a usable (finite, non-zero) H row
H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print('2. Run UINMF RPCA')
print('--------------------------------------')
data_rpca = integrate_rpca(data, H = H, W = W)

print('3. UINMF RPCA Postprocessing')
print('--------------------------------------')
data_rpca = data_rpca |>
  RunUMAP(uwot.init = 'pca', reduction = 'rpca_uinmf', dims = 1:min(20, round(dim(data_rpca[['rpca_uinmf']])[2]/2)),
          reduction.name = 'umap_rpca_uinmf')

print('4. UINMF RPCA Save data')
print('--------------------------------------')
SaveSeuratRds(data_rpca, str_glue('{out_path}04.16_RPCA_UINMF.rds'))
as.anndata(data_rpca, file_path = out_path, file_name = '04.16_RPCA_UINMF.h5ad',
           assay = 'Xenium', main_layer = 'data', other_layers = 'scale.data')

print('UINMF Seurat RPCA Done!')}
