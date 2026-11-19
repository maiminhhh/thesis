print('1. UINMF Harmony setup')
print('--------------------------------------')
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'fetal_liver'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
source('integration_datasets.R')

for (d in names(datasets)) {
print(str_glue('===== {d} ====='))
data = readRDS(str_glue('{datasets[[d]]$concat}.rds'))
out_path = datasets[[d]]$out
# integrate_harmony keeps only cells with a usable (finite, non-zero) H row
H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
W = readRDS(str_glue('{out_path}uinmf_W.rds'))
print('2. Run UINMF Harmony')
print('--------------------------------------')
data_harmony = integrate_harmony(data, H = H, W = W)
print('3. UINMF Harmony Postprocessing')
print('--------------------------------------')
data_harmony = data_harmony |>
  RunUMAP(uwot.init = 'pca', reduction = 'harmony',
          dims = 1:min(20, round(dim(data_harmony[['harmony']])[2]/2)),
          reduction.name = 'umap_harmony_uinmf')
print('3. UINMF Harmony Save data')
print('--------------------------------------')
SaveSeuratRds(data_harmony, str_glue('{out_path}04.14_Harmony_UINMF.rds'))
as.anndata(data_harmony, file_path = out_path,
           file_name = '04.14_Harmony_UINMF.h5ad',
           assay = 'Xenium', main_layer = 'counts', other_layers = c('data', 'scale.data'))
print('UINMF Harmony Done!')
}
