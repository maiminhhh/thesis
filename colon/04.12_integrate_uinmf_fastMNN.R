print('1. UINMF FastMNN Setup')
print('------------------------------------------------------')
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
source('integration_datasets.R')

for (d in names(datasets)) {
print(str_glue('===== {d} ====='))
data = readRDS(str_glue('{datasets[[d]]$concat}.rds'))
out_path = datasets[[d]]$out

# H comes from 04.11 and may cover fewer cells; integrate_fastmnn keeps only cells
# with a usable (finite, non-zero) H row
H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print('2. Run UINMF FastMNN')
print('------------------------------------------------------')
data_uinmf_fastmnn = integrate_fastmnn(data = data, H = H, W = W)

print('3. UINMF FastMNN Postprocessing')
print('------------------------------------------------------')
data_uinmf_fastmnn = data_uinmf_fastmnn |>
  RunUMAP(uwot.init = 'pca', reduction = 'mnn_uinmf', dims = 1:min(20, round(dim(data_uinmf_fastmnn[['mnn_uinmf']])[2]/2)),
          reduction.name = 'umap_mnn_uinmf')

print('4. UINMF FastMNN Save data')
print('------------------------------------------------------')
SaveSeuratRds(data_uinmf_fastmnn, str_glue('{out_path}04.12_FastMNN_UINMF.rds'))
as.anndata(data_uinmf_fastmnn, file_path = out_path,
           file_name = '04.12_FastMNN_UINMF.h5ad', assay = 'Xenium', main_layer = 'counts',
           other_layers = c('data', 'scale.data'))

print('UINMF FastMNN Done!')}
