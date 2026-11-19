print('1. UINMF CCA Set up')
print('--------------------------------------')
source('../library.R')
source('../source.R')
source('integration_sets.R')
param = MulticoreParam(workers = 32)

for (s in run_sets) {
concat = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print(str_glue('2. Run UINMF CCA - {s}'))
print('--------------------------------------')
data_cca = integrate_cca(data = concat, H = H, W = W)

print('3. UINMF CCA Postprocessing')
print('--------------------------------------')
# replace non-finite values in the integrated embedding before UMAP
red = data_cca[['cca_uinmf']]
emb = Embeddings(red)
emb[!is.finite(emb)] = 0

data_cca[['cca_uinmf']] = CreateDimReducObject(
  embeddings = emb,
  loadings   = Loadings(red),
  key        = Key(red),
  assay      = DefaultAssay(red)
)

data_cca = data_cca |>
  RunUMAP(reduction = 'cca_uinmf', dims = 1:min(20, round(dim(data_cca[['cca_uinmf']])[2]/2)),
          reduction.name = 'umap_cca_uinmf')

print('4. UINMF CCA Save data')
SaveSeuratRds(data_cca, str_glue('{out_path}04.13_CCA_UINMF.rds'))
as.anndata(data_cca, file_path = out_path, file_name = '04.13_CCA_UINMF.h5ad',
           assay = 'Xenium', main_layer = 'data', other_layers = 'scale.data')

print(str_glue('UINMF Seurat CCA Done! - {s}'))
}
