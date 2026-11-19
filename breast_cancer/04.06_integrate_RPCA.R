print('1. RPCA Set up')
print('--------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 32)
# data = readRDS('data/03_concatenated_data.rds')
# out_path = 'data/04_integrated_data/full_mouse/'

# print('2. Run RPCA')
# print('--------------------------------------')
# data_rpca = integrate_rpca(data)

# print('3. RPCA Postprocessing')
# print('--------------------------------------')
# DefaultAssay(data_rpca) = 'integrated'
# data_rpca = data_rpca |> 
#   ScaleData() |> 
#   RunPCA(reduction.name = 'rpca') |> 
#   RunUMAP(reduction = 'rpca', dims = 1:20,
#           reduction.name = 'umap_rpca')

# print('4. RPCA Save data')
# print('--------------------------------------')
# SaveSeuratRds(data_rpca, str_glue('{out_path}04.06_RPCA.rds'))
# as.anndata(data_rpca, file_path = out_path, file_name = '04.06_RPCA.h5ad',
#            assay = 'integrated', main_layer = 'data', other_layers = 'scale.data')

# print('Seurat RPCA Done!')

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.06_integrate_RPCA.R genes8` to run one dataset (datasets can run in parallel)
source('integration_sets.R')
source('umap_helper.R')   # RunUMAP, with init = 'pca' for genes6
for (i in run_sets) {

data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')

print('2. Run RPCA')
print('--------------------------------------')
data_rpca = integrate_rpca(data, half_dims = FALSE)

print('3. RPCA Postprocessing')
print('--------------------------------------')
# integrate_rpca already scales the integrated assay and builds the 'rpca' PCA (up to 20 PCs)
DefaultAssay(data_rpca) = 'integrated'
data_rpca = data_rpca |> 
  run_umap(i, 
  reduction = "rpca",
  dims = 1:min(20, ncol(data_rpca[['rpca']])),
  reduction.name = "umap_rpca",
  umap.method = "uwot",
  uwot.sgd = TRUE
)
print('4. RPCA Save data')
print('--------------------------------------')
SaveSeuratRds(data_rpca, str_glue('{out_path}04.06_RPCA.rds'))
write_h5ad_direct(data_rpca,
                  file  = str_glue('{out_path}04.06_RPCA.h5ad'),
                  assay = 'integrated',
                  layer = 'data',
                  reductions = c('rpca', 'umap_rpca'))
print('Seurat RPCA Done!')
}