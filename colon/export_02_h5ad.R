# Check that the 02 objects carry the EBCC labels and write their h5ad files
# (as.anndata needs RETICULATE_PYTHON pointing at an env with anndata)
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
for (s in c('human', 'mouse')) {
  x = readRDS(sprintf('data/02_annotated_data/02_%s_annotated.rds', s))
  stopifnot('ebcc_label' %in% colnames(x[[]]), !anyNA(x$ebcc_label))
  print(sprintf('%s: %d cells, EBCC labels:', s, ncol(x))); print(table(x$ebcc_label))
  as.anndata(x, file_path = 'data/02_annotated_data/',
             file_name = sprintf('02_%s_annotated.h5ad', s), assay = 'Xenium',
             main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')
}
