print('1. FastMNN Set up')
print('--------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 32)
# data = readRDS('data/03_concatenated_data.rds')
# out_path = 'data/04_integrated_data/full_mouse/'

# print('2. Run FastMNN')
# print('--------------------------------------')
# data_fastmnn = integrate_fastmnn(data)
# data_fastmnn = JoinLayers(data_fastmnn)
# data_fastmnn = data_fastmnn |> 
#   RunUMAP(reduction = 'mnn', dims = 1:20,
#           reduction.name = 'umap_mnn')

# print('3. FastMNN Save data')
# print('--------------------------------------')
# SaveSeuratRds(data_fastmnn, str_glue('{out_path}04.01_FastMNN.rds'))

# as.anndata(data_fastmnn, file_path = out_path, file_name = '04.01_FastMNN.h5ad', 
#            assay = 'MERFISH', 
#            main_layer = 'counts', other_layers = c('data', 'scale.data'))

# print('FastMNN Done!')

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.01_integrate_fastMNN.R genes8` to run one dataset (datasets can run in parallel)
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {
data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')

print('2. Run FastMNN')
print('--------------------------------------')
data_fastmnn = integrate_fastmnn(data, half_dims = FALSE)
data_fastmnn = JoinLayers(data_fastmnn)
data_fastmnn = data_fastmnn |> 
  run_umap(i, reduction = 'mnn', dims = 1:min(20, ncol(data_fastmnn[['mnn']])),
          reduction.name = 'umap_mnn')

print('3. FastMNN Save data')
print('--------------------------------------')
SaveSeuratRds(data_fastmnn, str_glue('{out_path}04.01_FastMNN.rds'))

as.anndata(data_fastmnn, file_path = out_path, file_name = '04.01_FastMNN.h5ad', 
           assay = 'MERFISH', 
           main_layer = 'counts', other_layers = c('data', 'scale.data'))

print('FastMNN Done!')
}