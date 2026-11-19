# attach the scVI label transfer (06_label_transfer_scVI.py, human -> mouse, dissected)
# to the original-panel objects, so downstream analysis can use all panel genes.
#   mouse (dissected): celltype_transfer = label transferred from human,
#                      transfer_confidence = share of the 30 human neighbours voting for it
#   human:             celltype_transfer = ebcc_label (the reference labels used)
# Cells not in the scVI run (no counts on the shared orthologues) get NA.
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
out_dir = 'data/06_label_transfer/dissected/'

transfer = read.csv(str_glue('{out_dir}06_scVI_human_to_mouse_labels.csv'), row.names = 1)
rownames(transfer) = sub('^mmusculus_', '', rownames(transfer))
ref_cells = readLines(str_glue('{out_dir}06_scVI_human_to_mouse_reference_cells.txt'))
print('Done setting up!')

# mouse (dissected)
mouse = readRDS('data/02_annotated_data/02_mouse_subset_annotated.rds')
mouse$celltype_transfer = transfer[colnames(mouse), 'after']
mouse$transfer_confidence = transfer[colnames(mouse), 'confidence']
print(str_glue('mouse: {sum(!is.na(mouse$celltype_transfer))} / {ncol(mouse)} cells labelled'))
print(table(ebcc = mouse$ebcc_label, transfer = mouse$celltype_transfer, useNA = 'ifany'))

# human: reference labels, only for the cells that were in the scVI run
human = readRDS('data/02_annotated_data/02_human_annotated.rds')
in_scvi = colnames(human) %in% sub('^hsapiens_', '', ref_cells)
human$celltype_transfer = ifelse(in_scvi, as.character(human$ebcc_label), NA)
print(str_glue('human: {sum(in_scvi)} / {ncol(human)} cells in the scVI run'))

for (x in list(list(mouse, '06_mouse_dissected_transfer'), list(human, '06_human_transfer'))) {
  SaveSeuratRds(x[[1]], str_glue('{out_dir}{x[[2]]}.rds'))
  as.anndata(x[[1]], file_path = out_dir, file_name = x[[2]], assay = 'Xenium',
             main_layer = 'counts', other_layers = c('data', 'scale.data'))
  # cell centroids + labels for the spatial analysis (07); the h5ad carries no coordinates
  xy = GetTissueCoordinates(x[[1]])
  meta = x[[1]]@meta.data[xy$cell, intersect(c('ebcc_label', 'celltype_transfer', 'transfer_confidence'),
                                             colnames(x[[1]]@meta.data))]
  write.csv(cbind(xy[, c('cell', 'x', 'y')], meta), str_glue('{out_dir}{x[[2]]}_cells.csv'), row.names = FALSE)
}

plot_mouse_spatial = function(mouse) {
  # same colours as the annotation spatial plots (hue palette over the sorted labels)
  labs = sort(unique(na.omit(c(as.character(mouse$ebcc_label), mouse$celltype_transfer))))
  cols = setNames(scales::hue_pal()(length(labs)), labs)
  # same factor levels in both columns, so the two panels share one legend
  for (col in c('ebcc_label', 'celltype_transfer'))
    mouse@meta.data[[col]] = factor(as.character(mouse@meta.data[[col]]), levels = labs)
  titles = c(ebcc_label = 'mouse EBCC label (before)',
             celltype_transfer = 'transferred from human, scVI (after)')
  plots = Map(\(col, title) ImageDimPlot(mouse, fov = Images(mouse)[1], group.by = col, cols = cols,
                                         size = 0.8, axes = TRUE) + labs(title = title),
              names(titles), titles)
  p = patchwork::wrap_plots(plots, ncol = 2, guides = 'collect') & labs(fill = 'cell type')
  dir.create('../figures/bone_marrow/label_transfer', showWarnings = FALSE, recursive = TRUE)
  ggsave('../figures/bone_marrow/label_transfer/spatial_mouse_dissected.png', p,
         width = 12, height = 12, dpi = 200, bg = 'white')
}
plot_mouse_spatial(mouse)
print('Done!')
