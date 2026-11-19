print('1. scVI-latent FastMNN Set up')
print('------------------------------------------------------')
# FastMNN on the scVI latent (04.18, no species key) instead of PCA
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
source('integration_sets.R')

for (s in run_sets) {
data = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
Z = read_scvi_latent(out_path)

print(str_glue('2. Run scVI-latent FastMNN - {s}'))
print('------------------------------------------------------')
data_fastmnn = integrate_fastmnn(data = data, H = Z, W = matrix(numeric(0), 0, 0)) |>
  rename_reduction('mnn_uinmf', 'mnn_scvi')

print('3. scVI-latent FastMNN Postprocessing')
print('------------------------------------------------------')
data_fastmnn = data_fastmnn |>
  RunUMAP(reduction = 'mnn_scvi', dims = 1:min(20, round(dim(data_fastmnn[['mnn_scvi']])[2]/2)),
          reduction.name = 'umap_mnn_scvi')

print('4. scVI-latent FastMNN Save data')
print('------------------------------------------------------')
SaveSeuratRds(data_fastmnn, str_glue('{out_path}04.19_FastMNN_scVI.rds'))
as.anndata(data_fastmnn, file_path = out_path,
           file_name = '04.19_FastMNN_scVI.h5ad', assay = 'Xenium', main_layer = 'counts',
           other_layers = c('data', 'scale.data'))

print(str_glue('scVI-latent FastMNN Done! - {s}'))
}
