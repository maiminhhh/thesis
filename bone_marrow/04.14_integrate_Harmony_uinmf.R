print('1. UINMF Harmony setup')
print('--------------------------------------')
source('../library.R')
source('../source.R')
source('integration_sets.R')
param = MulticoreParam(workers = 32)

for (s in run_sets) {
data = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print(str_glue('2. Run UINMF Harmony - {s}'))
print('--------------------------------------')
data_harmony = integrate_harmony(data, H = H, W = W)

print('3. UINMF Harmony Postprocessing')
print('--------------------------------------')
data_harmony = data_harmony |>
  RunUMAP(reduction = 'harmony',
          dims = 1:min(20, ncol(data_harmony[['harmony']])),
          reduction.name = 'umap_harmony_uinmf')

print('4. UINMF Harmony Save data')
print('--------------------------------------')
SaveSeuratRds(data_harmony, str_glue('{out_path}04.14_Harmony_UINMF.rds'))
as.anndata(data_harmony, file_path = out_path,
           file_name = '04.14_Harmony_UINMF.h5ad',
           assay = 'Xenium', main_layer = 'counts', other_layers = c('data', 'scale.data'))

print(str_glue('UINMF Harmony Done! - {s}'))
}
