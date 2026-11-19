print('1. UINMF RPCA Set up')
print('--------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
source('integration_sets.R')
param = MulticoreParam(workers = 32)

for (s in run_sets) {
data = readRDS(integration_sets[[s]]$concat)
out_path = integration_sets[[s]]$out
H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))

print(str_glue('2. Run UINMF RPCA - {s}'))
print('--------------------------------------')
data_rpca = integrate_rpca(data, H = H, W = W)

print('3. UINMF RPCA Postprocessing')
print('--------------------------------------')
d = min(20, round(dim(data_rpca[['rpca_uinmf']])[2]/2))
data_rpca = data_rpca |>
  RunUMAP(reduction = 'rpca_uinmf', dims = 1:d,
          reduction.name = 'umap_rpca_uinmf')

print('4. UINMF RPCA Save data')
print('--------------------------------------')
SaveSeuratRds(data_rpca, str_glue('{out_path}04.16_RPCA_UINMF.rds'))
write_h5ad_direct(data_rpca,
                  file  = str_glue('{out_path}04.16_RPCA_UINMF.h5ad'),
                  assay = DefaultAssay(data_rpca),
                  layer = 'data',
                  reductions = Reductions(data_rpca))

print(str_glue('UINMF Seurat RPCA Done! - {s}'))
}
