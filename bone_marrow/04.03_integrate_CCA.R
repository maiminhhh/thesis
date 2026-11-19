print('1. CCA Set up')
print('--------------------------------------')
source('../library.R')
source('../source.R')
source('integration_sets.R')
param = MulticoreParam(workers = 32)
options(future.globals.maxSize = 20 * 1024^3)

for (s in run_sets) {
print(str_glue('2. Run CCA - {s}'))
concat = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
data_cca = integrate_cca(data = concat)
data_cca = data_cca |>
  ScaleData() |>
  RunPCA(reduction.name = 'cca') |>
  RunUMAP(reduction = 'cca', dims = 1:min(20, round(dim(data_cca[['cca']])[2]/2)), reduction.name = 'umap_cca')

print('3. CCA Save data')
SaveSeuratRds(data_cca, str_glue('{out_path}04.03_CCA.rds'))
as.anndata(data_cca, file_path = out_path, file_name = '04.03_CCA.h5ad',
           assay = 'integrated', main_layer = 'data', other_layers = 'scale.data')
}
