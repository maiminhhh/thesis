
print('1. UINMF RPCA Set up')
print('--------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 32)
# data = readRDS('data/03_concatenated_data.rds')
# out_path = 'data/04_integrated_data/full_mouse/'

# H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
# W = readRDS(str_glue('{out_path}uinmf_W.rds'))

# print('2. Run UINMF RPCA')
# print('--------------------------------------')
# data_rpca = integrate_rpca(data, H = H, W = W)

# print('3. UINMF RPCA Postprocessing')
# print('--------------------------------------')
# data_rpca = data_rpca |> 
#   RunUMAP(reduction = 'rpca_uinmf', dims = 1:20,
#           reduction.name = 'umap_rpca_uinmf')

# print('4. UINMF RPCA Save data')
# print('--------------------------------------')
# SaveSeuratRds(data_rpca, str_glue('{out_path}04.16_RPCA_UINMF.rds'))
# #write_h5ad(data_rpca, str_glue('{out_path}04.16_RPCA_UINMF.h5ad'))
# as.anndata(data_rpca, file_path = out_path, file_name = '04.16_RPCA_UINMF.h5ad',
#            assay = 'MERFISH', main_layer = 'data', other_layers = 'scale.data')

# print('UINMF Seurat RPCA Done!')

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.16_integrate_uinmf_RPCA.R genes8` to run one dataset (datasets can run in parallel)
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {
data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')

H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print('2. Run UINMF RPCA')
print('--------------------------------------')
data_rpca = integrate_rpca(data, H = H, W = W, half_dims = FALSE)

print('3. UINMF RPCA Postprocessing')
print('--------------------------------------')
data_rpca = data_rpca |> 
  run_umap(i, reduction = 'rpca_uinmf', dims = 1:min(20, ncol(data_rpca[['rpca_uinmf']])),
          reduction.name = 'umap_rpca_uinmf')

print('4. UINMF RPCA Save data')
print('--------------------------------------')
SaveSeuratRds(data_rpca, str_glue('{out_path}04.16_RPCA_UINMF.rds'))
#write_h5ad(data_rpca, str_glue('{out_path}04.16_RPCA_UINMF.h5ad'))
as.anndata(data_rpca, file_path = out_path, file_name = '04.16_RPCA_UINMF.h5ad',
           assay = 'MERFISH', main_layer = 'data', other_layers = 'scale.data')

print('UINMF Seurat RPCA Done!')}