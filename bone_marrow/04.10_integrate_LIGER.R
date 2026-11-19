print('1. LIGER Set up')
print('--------------------------------------')
source('../library.R')
source('../source.R')
source('integration_sets.R')
param = MulticoreParam(workers = 32)

for (s in run_sets) {
data = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
print(str_glue('2. Run LIGER - {s}'))
data_liger = integrate_liger(data = data)

print('3. LIGER Postprocessing')
data_liger = data_liger |>
  RunUMAP(reduction = 'inmfNorm', dims = 1:20,
          reduction.name = 'umap_liger')

print('4. LIGER Save data')
print('--------------------------------------')
SaveSeuratRds(data_liger, str_glue('{out_path}04.10_LIGER.rds'))
as.anndata(data_liger, file_path = out_path, file_name = '04.10_LIGER.h5ad',
           assay = 'Xenium', main_layer = 'counts', other_layers = c('data', 'scale.data'))
print(str_glue('LIGER Done! - {s}'))
}
