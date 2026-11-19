print('1. scVI-latent RPCA Set up')
print('--------------------------------------')
# Seurat RPCA anchors on the scVI latent (04.18, no species key) instead of PCA
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'fetal_liver'))
source('../library.R')
source('../source.R')
source('integration_datasets.R')

for (s in names(datasets)) {
print(str_glue("===== {s} ====="))
data = readRDS(str_glue("{datasets[[s]]$concat}.rds"))
out_path = datasets[[s]]$out
Z = read_scvi_latent(out_path)
# IntegrateEmbeddings carries the input reduction's gene loadings along (it does not
# use them), and needs one column per dimension: the latent has none, so zeros
W = matrix(0, nrow(data), ncol(Z$hsapiens), dimnames = list(rownames(data), NULL))

print(str_glue('2. Run scVI-latent RPCA - {s}'))
print('--------------------------------------')
data_rpca = integrate_rpca(data, H = Z, W = W) |>
  rename_reduction('rpca_uinmf', 'rpca_scvi') |>
  rename_reduction('uinmf_H', 'scvi_latent')

print('3. scVI-latent RPCA Postprocessing')
print('--------------------------------------')
d = min(20, round(dim(data_rpca[['rpca_scvi']])[2]/2))
data_rpca = data_rpca |>
  RunUMAP(uwot.init = 'pca', reduction = 'rpca_scvi', dims = 1:d,
          reduction.name = 'umap_rpca_scvi')

print('4. scVI-latent RPCA Save data')
print('--------------------------------------')
SaveSeuratRds(data_rpca, str_glue('{out_path}04.23_RPCA_scVI.rds'))
write_h5ad_direct(data_rpca,
                  file  = str_glue('{out_path}04.23_RPCA_scVI.h5ad'),
                  assay = DefaultAssay(data_rpca),
                  layer = 'data',
                  reductions = Reductions(data_rpca))

print(str_glue('scVI-latent Seurat RPCA Done! - {s}'))
}
