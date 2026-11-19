setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'bone_marrow'))

print('1. Harmony setup')
print('--------------------------------------')
source('../library.R')
source('../source.R')
source('integration_sets.R')
param = MulticoreParam(workers = 32)

for (s in run_sets) {
data = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
print(str_glue('2. Run Harmony - {s}'))
print('--------------------------------------')
data_harmony = integrate_harmony(data)
# harmony has as many dimensions as PCs, which is capped by the number of genes
data_harmony = data_harmony |>
  RunUMAP(reduction = 'harmony', dims = 1:min(20, ncol(data_harmony[['harmony']])),
          reduction.name = 'umap_harmony')

print('3. Save data')
print('--------------------------------------')
SaveSeuratRds(data_harmony, str_glue('{out_path}04.04_Harmony.rds'))
as.anndata(data_harmony, file_path = out_path, file_name = '04.04_Harmony.h5ad',
           assay = 'Xenium', main_layer = 'counts', other_layers = c('data', 'scale.data'))

print(str_glue('Harmony Done! - {s}'))
}
