
print('1. UINMF CCA Set up')
print('--------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 32)
# concat = readRDS('data/03_concatenated_data.rds')
# out_path = 'data/04_integrated_data/full_mouse/'

# H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
# W = readRDS(str_glue('{out_path}uinmf_W.rds'))

# print('3. Run UINMF CCA')
# print('--------------------------------------')
# data_cca = integrate_cca(data = concat, H = H, W = W)

# print('4. UINMF CCA Postprocessing')
# print('--------------------------------------')
# data_cca = data_cca |> 
#   RunUMAP(reduction = 'cca_uinmf', dims = 1:20,
#           reduction.name = 'umap_cca_uinmf')

# print('4. UINMF CCA Save data')

# SaveSeuratRds(data_cca, str_glue('{out_path}04.13_CCA_UINMF.rds'))
# #write_h5ad(data_cca, str_glue('{out_path}04.13_CCA_UINMF.h5ad'))
# as.anndata(data_cca, file_path = out_path, file_name = '04.13_UINMF_CCA.h5ad', 
#            assay = 'MERFISH', main_layer = 'data', other_layers = 'scale.data')

# print('UINMF Seurat CCA Done!')

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.13_integrate_uinmf_CCA.R genes8` to run one dataset (datasets can run in parallel)
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {
concat = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')

H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print('3. Run UINMF CCA')
print('--------------------------------------')
data_cca = integrate_cca(data = concat, H = H, W = W, half_dims = FALSE)

print('4. UINMF CCA Postprocessing')
print('--------------------------------------')
data_cca = data_cca |> 
  run_umap(i, reduction = 'cca_uinmf', dims = 1:min(20, ncol(data_cca[['cca_uinmf']])),
          reduction.name = 'umap_cca_uinmf')

print('4. UINMF CCA Save data')

SaveSeuratRds(data_cca, str_glue('{out_path}04.13_CCA_UINMF.rds'))
#write_h5ad(data_cca, str_glue('{out_path}04.13_CCA_UINMF.h5ad'))
as.anndata(data_cca, file_path = out_path, file_name = '04.13_CCA_UINMF.h5ad', 
           assay = 'MERFISH', main_layer = 'data', other_layers = 'scale.data')

print('UINMF Seurat CCA Done!')}