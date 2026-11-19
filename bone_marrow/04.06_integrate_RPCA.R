print('1. RPCA Set up')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
print('--------------------------------------')
source('../library.R')
source('../source.R')
source('integration_sets.R')
param = MulticoreParam(workers = 32)

for (s in run_sets) {
print(str_glue('2. Run RPCA - {s}'))
data = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
data_rpca = integrate_rpca(data)

print('3. RPCA Postprocessing')
print('--------------------------------------')
DefaultAssay(data_rpca) = 'integrated'
d = min(20, round(dim(data_rpca[['rpca']])[2]/2))
# UMAP_INIT=random avoids a very slow spectral start when many cells share
# near-identical embeddings (genes2); the UMAP is only used for plots
data_rpca = data_rpca |>
  RunUMAP(reduction = 'rpca', dims = 1:d,
          reduction.name = 'umap_rpca',
          uwot.init = Sys.getenv('UMAP_INIT', 'spectral'))

print('4. RPCA Save data')
print('--------------------------------------')
SaveSeuratRds(data_rpca, str_glue('{out_path}04.06_RPCA.rds'))
# the integrated assay has no counts layer (as.anndata's default other_layers)
as.anndata(data_rpca, file_path = out_path, file_name = '04.06_RPCA.h5ad',
           assay = 'integrated', main_layer = 'data', other_layers = 'scale.data')

print(str_glue('Seurat RPCA Done! - {s}'))
}
