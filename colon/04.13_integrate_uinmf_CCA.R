print('1. UINMF CCA Set up')
print('--------------------------------------')
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
source('integration_datasets.R')

for (d in names(datasets)) {
print(str_glue('===== {d} ====='))
concat = readRDS(str_glue('{datasets[[d]]$concat}.rds'))
out_path = datasets[[d]]$out

# H comes from 04.11 and may cover fewer cells; integrate_cca keeps only cells
# with a usable (finite, non-zero) H row
H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print('3. Run UINMF CCA')
print('--------------------------------------')
data_cca = integrate_cca(data = concat, H = H, W = W)

print('4. UINMF CCA Postprocessing')
print('--------------------------------------')
# the integrated embedding can carry a few non-finite values, which break UMAP and scib
red = data_cca[['cca_uinmf']]
emb = Embeddings(red)
print(sprintf('non-finite values in cca_uinmf set to 0: %d', sum(!is.finite(emb))))
emb[!is.finite(emb)] = 0
data_cca[['cca_uinmf']] = CreateDimReducObject(embeddings = emb, loadings = Loadings(red),
                                               key = Key(red), assay = DefaultAssay(red))
data_cca = data_cca |>
  RunUMAP(uwot.init = 'pca', reduction = 'cca_uinmf', dims = 1:min(20, round(dim(data_cca[['cca_uinmf']])[2]/2)),
          reduction.name = 'umap_cca_uinmf')

print('4. UINMF CCA Save data')
SaveSeuratRds(data_cca, str_glue('{out_path}04.13_CCA_UINMF.rds'))
# as.anndata fails on the split-then-joined assay this leaves: write the h5ad directly
write_h5ad_direct(data_cca, file = str_glue('{out_path}04.13_CCA_UINMF.h5ad'),
                  assay = 'Xenium', layer = 'data', reductions = Reductions(data_cca))

print('UINMF Seurat CCA Done!')}
