print('1. UINMF FastMNN Setup')
print('------------------------------------------------------')
source('../library.R')
source('../source.R')
source('integration_sets.R')
param = MulticoreParam(workers = 32)

for (s in run_sets) {
data = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print(str_glue('2. Run UINMF FastMNN - {s}'))
print('------------------------------------------------------')
data_uinmf_fastmnn = integrate_fastmnn(data = data, H = H, W = W)

print('3. UINMF FastMNN Postprocessing')
print('------------------------------------------------------')
data_uinmf_fastmnn = data_uinmf_fastmnn |>
  RunUMAP(reduction = 'mnn_uinmf', dims = 1:min(20, round(dim(data_uinmf_fastmnn[['mnn_uinmf']])[2]/2)),
          reduction.name = 'umap_mnn_uinmf')

print('4. UINMF FastMNN Save data')
print('------------------------------------------------------')
SaveSeuratRds(data_uinmf_fastmnn, str_glue('{out_path}04.12_FastMNN_UINMF.rds'))
as.anndata(data_uinmf_fastmnn, file_path = out_path,
           file_name = '04.12_FastMNN_UINMF.h5ad', assay = 'Xenium', main_layer = 'counts',
           other_layers = c('data', 'scale.data'))

print(str_glue('UINMF FastMNN Done! - {s}'))
}
