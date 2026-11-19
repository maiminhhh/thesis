print('1. FastMNN Set up')
print('--------------------------------------')
source('../library.R')
source('../source.R')
source('integration_sets.R')
param = MulticoreParam(workers = 32)

for (s in run_sets) {
data = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
print(str_glue('2. Run FastMNN - {s}'))
print('--------------------------------------')
data_fastmnn = integrate_fastmnn(data)
data_fastmnn = JoinLayers(data_fastmnn)

# mnn has at most 10 dimensions (see integrate_fastmnn)
data_fastmnn = data_fastmnn |>
  RunUMAP(reduction = 'mnn', dims = 1:ncol(data_fastmnn[['mnn']]),
          reduction.name = 'umap_mnn')

print('3. FastMNN Save data')
print('--------------------------------------')
SaveSeuratRds(data_fastmnn, str_glue('{out_path}04.01_FastMNN.rds'))
as.anndata(data_fastmnn, file_path = out_path, file_name = '04.01_FastMNN.h5ad',
           assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))

print(str_glue('FastMNN Done! - {s}'))
}
