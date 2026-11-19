# genes6 variant of 04.03_integrate_CCA.R (figure-only UMAP with init = 'pca')
print('1. CCA Set up')
print('--------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 32)
# concat = readRDS('data/03_concatenated_data.rds')
# out_path = 'data/04_integrated_data/full_mouse/'
# options(future.globals.maxSize = 20 * 1024^3)

# print('2. Run CCA')
# print('--------------------------------------')
# data_cca = integrate_cca(data = concat)
# data_cca = data_cca |> 
#   ScaleData() |> 
#   RunPCA(reduction.name = 'cca') |> 
#   RunUMAP(reduction = 'cca', dims = 1:20, reduction.name = 'umap_cca')

# print('3. CCA Save data')
# SaveSeuratRds(data_cca, str_glue('{out_path}04.03_CCA.rds'))
# write_h5ad(data_cca, str_glue('{out_path}04.03_CCA.h5ad'))
# as.anndata(data_cca, file_path = out_path, file_name = '04.03_CCA.h5ad', 
#            assay = 'integrated', main_layer = 'data', other_layers = 'scale.data')

# print('Seurat CCA Done!')

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.03_integrate_CCA.R genes8` to run one dataset (datasets can run in parallel)
source('integration_sets.R')

# genes6 only: RunUMAP's spectral initialisation crashes on this dataset (disconnected kNN graph with 49 genes),
# so the figure-only UMAP uses Seurat's RunUMAP settings with init = 'pca'; the integration is unchanged
umap_pca_init = function(obj, reduction, dims, name) {
  set.seed(42)
  emb = uwot::umap(Embeddings(obj, reduction)[, dims, drop = FALSE], n_neighbors = 30, metric = 'cosine',
                   min_dist = 0.3, init = 'pca', n_threads = as.integer(Sys.getenv('N_THREADS', 4)))
  rownames(emb) = colnames(obj); colnames(emb) = paste0(gsub('_', '', name), '_', seq_len(ncol(emb)))
  obj[[name]] = CreateDimReducObject(embeddings = emb, key = paste0(gsub('_', '', name), '_'), assay = DefaultAssay(obj))
  obj
}
for (i in run_sets) {
concat = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))
out_path = str_glue('data/{i}/')
options(future.globals.maxSize = 20 * 1024^3)

print('2. Run CCA')
print('--------------------------------------')
# integrate_cca already scales the integrated assay and builds the 'cca' PCA (up to 20 PCs)
data_cca = integrate_cca(data = concat)
data_cca = data_cca |> 
  umap_pca_init('cca', 1:min(20, ncol(data_cca[['cca']])), 'umap_cca')

print('3. CCA Save data')
SaveSeuratRds(data_cca, str_glue('{out_path}04.03_CCA.rds'))
# write_h5ad(data_cca, str_glue('{out_path}04.03_CCA.h5ad'))
as.anndata(data_cca, file_path = out_path, file_name = '04.03_CCA.h5ad', 
           assay = 'integrated', main_layer = 'data', other_layers = 'scale.data')

print('Seurat CCA Done!')
}