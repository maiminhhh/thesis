print('1. UINMF FastMNN Setup')
print('------------------------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 32)
# data = readRDS('data/03_concatenated_data.rds')
# out_path = 'data/04_integrated_data/full_mouse/'

# H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
# W = readRDS(str_glue('{out_path}uinmf_W.rds'))

# print('2. Run UINMF FastMNN')
# print('------------------------------------------------------')
# data_uinmf_fastmnn = integrate_fastmnn(data = data, H = H, W = W)

# print('3. UINMF FastMNN Postprocessing')
# print('------------------------------------------------------')
# data_uinmf_fastmnn = data_uinmf_fastmnn |> 
#   RunUMAP(reduction = 'mnn_uinmf', dims = 1:20,
#           reduction.name = 'umap_mnn_uinmf')

# print('4. UINMF FastMNN Save data')
# print('------------------------------------------------------')
# SaveSeuratRds(data_uinmf_fastmnn, str_glue('{out_path}04.12_FastMNN_UINMF.rds'))
# #write_h5ad(data_uinmf_fastmnn, str_glue('{out_path}04.12_FastMNN_UINMF.h5ad'))
# as.anndata(data_uinmf_fastmnn, file_path = out_path, 
#            file_name = '04.12_FastMNN_UINMF.h5ad', assay = 'MERFISH', main_layer = 'counts',
#            other_layers = c('data', 'scale.data'))

# print('UINMF FastMNN Done!')

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.12_integrate_uinmf_fastMNN.R genes8` to run one dataset (datasets can run in parallel)
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {
data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')

H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print('2. Run UINMF FastMNN')
print('------------------------------------------------------')
data_uinmf_fastmnn = integrate_fastmnn(data = data, H = H, W = W)

print('3. UINMF FastMNN Postprocessing')
print('------------------------------------------------------')
data_uinmf_fastmnn = data_uinmf_fastmnn |> 
  run_umap(i, reduction = 'mnn_uinmf', dims = 1:min(20, ncol(data_uinmf_fastmnn[['mnn_uinmf']])),
          reduction.name = 'umap_mnn_uinmf')

print('4. UINMF FastMNN Save data')
print('------------------------------------------------------')
SaveSeuratRds(data_uinmf_fastmnn, str_glue('{out_path}04.12_FastMNN_UINMF.rds'))
#write_h5ad(data_uinmf_fastmnn, str_glue('{out_path}04.12_FastMNN_UINMF.h5ad'))
as.anndata(data_uinmf_fastmnn, file_path = out_path, 
           file_name = '04.12_FastMNN_UINMF.h5ad', assay = 'MERFISH', main_layer = 'counts',
           other_layers = c('data', 'scale.data'))

print('UINMF FastMNN Done!')}