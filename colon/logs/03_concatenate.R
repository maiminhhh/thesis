setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
human = readRDS('data/02_annotated_data/02_human_annotated.rds')
homology_tbl = read.csv('../homology_tbl.csv')

# ImageDimPlot(mouse, group.by = 'celltype', axes = TRUE)
# ImageDimPlot(human, group.by = 'celltype', axes = TRUE)

# coords already has a 'cell' column; its rownames are just 1..n
coords = GetTissueCoordinates(mouse)
coords$celltype = mouse$celltype[match(coords$cell, colnames(mouse))]
# dissection (../EBCC_PIPELINE.md 5b): the human sample is a full-thickness non-diseased colon section, and the
# whole mouse section is full-thickness colon wall (no squamous, skeletal muscle, cartilage/bone; ~1% fat), so
# this rectangle is a same-tissue subset; options compared in logs/dissection/ (kept 2026-10-24)
xmin = 260;  xmax = 2685
ymin = 2590; ymax = 4625

# p = plot_ly(
#   coords,
#   x = ~x,
#   y = ~y,
#   color = ~celltype,
#   type = "scattergl",
#   mode = "markers",
#   marker = list(size = 3),
#   text = ~paste0(
#     "x: ", x,
#     "<br>y: ", y
#   ),
#   hovertemplate = "%{text}<extra></extra>"
# ) |>
#   layout(
#     dragmode = "drawclosedpath",
#     newshape = list(line = list(width = 2))
#   ) |>
#   config(
#     modeBarButtonsToAdd = c(
#       "drawclosedpath",
#       "eraseshape"
#     )
#   )
# 
# p
# 
# ggplot(coords, aes(x, y, color = celltype)) +
#   geom_point(size = 0.025) +
#   annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
#            fill = NA, color = "black", linewidth = 0.5) +
#   coord_fixed() +
#   theme_void() +
#   guides(color = guide_legend(override.aes = list(size = 3)))

coords$keep = coords$x >= xmin &
  coords$x <= xmax &
  coords$y >= ymin &
  coords$y <= ymax

mouse_subset = mouse |>
  subset(cells = coords$cell[coords$keep])
SaveSeuratRds(mouse_subset, 'data/02_annotated_data/02_mouse_subset_annotated.rds')
as.anndata(mouse_subset, file_path = 'data/02_annotated_data/',
           file_name = '02_mouse_subset_annotated.h5ad', assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))
print(sprintf('Dissected mouse: %d of %d cells', ncol(mouse_subset), ncol(mouse)))

# ImageDimPlot(mouse_subset, group.by = 'celltype', axes = TRUE)
# ImageDimPlot(human, group.by = 'celltype', axes = TRUE)

# Ensembl ids: concatenate() maps the species through the homology table on them
rownames(human) = human[['Xenium']]@meta.data$ensembl_id
rownames(mouse) = mouse[['Xenium']]@meta.data$ensembl_id
rownames(mouse_subset) = mouse_subset[['Xenium']]@meta.data$ensembl_id

# per-species inputs for SAMap / LIGER UINMF (Ensembl ids), restricted to the cells
# kept in the concatenated data, saved as 02_{human,mouse}_annotated in out_dir
save_species_inputs = function(concat, human, mouse, out_dir) {
  # concat cell names are '<species>_<cell id>' (add.cell.ids in concatenate)
  kept = colnames(concat)
  human_cells = sub('^hsapiens_', '', kept[concat$orig.ident == 'hsapiens'])
  mouse_cells = sub('^mmusculus_', '', kept[concat$orig.ident == 'mmusculus'])
  stopifnot(all(human_cells %in% colnames(human)), all(mouse_cells %in% colnames(mouse)))
  human = human[, human_cells]
  mouse = mouse[, mouse_cells]
  stopifnot(ncol(human) + ncol(mouse) == ncol(concat))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  SaveSeuratRds(human, file.path(out_dir, '02_human_annotated.rds'))
  SaveSeuratRds(mouse, file.path(out_dir, '02_mouse_annotated.rds'))
  as.anndata(human, file_path = out_dir,
             file_name = '02_human_annotated.h5ad', assay = 'Xenium',
             main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')
  as.anndata(mouse, file_path = out_dir,
             file_name = '02_mouse_annotated.h5ad', assay = 'Xenium',
             main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')
}

# drop cells with no counts on the shared genes, then scale both species jointly
# (merge() keeps each species' own scale.data) so the baseline stays unintegrated
prepare_concat = function(concat) {
  VariableFeatures(concat) = rownames(concat)
  concat = JoinLayers(concat)
  concat$nCount_Xenium = colSums(x = concat[['Xenium']]$counts)
  concat = concat |>
    subset(nCount_Xenium > 0)
  concat |>
    ScaleData(features = rownames(concat)) |>
    RunPCA(features = rownames(concat), npcs = min(50, nrow(concat) - 1)) |>
    # PCA init: with few genes many cells are identical, the kNN graph falls apart and
    # uwot's spectral init can hang (the UMAP here is only for viewing)
    RunUMAP(dims = 1:min(20, nrow(concat) - 1), uwot.init = 'pca')
}

save_concat = function(concat, file) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  SaveSeuratRds(concat, paste0(file, '.rds'))
  as.anndata(concat, file_path = paste0(dirname(file), '/'),
             file_name = paste0(basename(file), '.h5ad'), assay = 'Xenium',
             main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')
}

# human + full mouse
concat = prepare_concat(concatenate(human, mouse, homology_tbl))
save_concat(concat, 'data/03_concatenated_data/03_concatenated_data_full')
save_species_inputs(concat, human, mouse, 'data/03_concatenated_data/full_mouse/')
print(sprintf('full: %d genes x %d cells', nrow(concat), ncol(concat)))

# human + dissected mouse
concat = prepare_concat(concatenate(human, mouse_subset, homology_tbl))
save_concat(concat, 'data/03_concatenated_data/03_concatenated_data_subset')
save_species_inputs(concat, human, mouse_subset, 'data/03_concatenated_data/dissected_mouse/')
print(sprintf('dissected: %d genes x %d cells', nrow(concat), ncol(concat)))

# from human + dissected mouse; p = fraction of shared genes DROPPED, saved to data/genes<10p>/
shared_genes = concat[['Xenium']][[]]
# mouse id for every shared (human id) gene, looked up by name rather than by row order
hs_to_mm = setNames(shared_genes$mmusculus_homolog_ensembl_gene,
                    shared_genes$hsapiens_homolog_ensembl_gene)
stopifnot(all(rownames(concat) %in% names(hs_to_mm)))
n_shared = nrow(concat)
# nested drops: one random gene order, each level drops a longer prefix of it
set.seed(1)
gene_order = sample(rownames(concat))
# levels run one after another: forking them (mclapply) segfaults in RunPCA / RunUMAP
for (frac in c(0.2, 0.4, 0.6, 0.8)) {
  drop = gene_order[seq_len(round(frac * n_shared))]
  # each species renormalised (and rescaled, as in 02) on its remaining genes, then concatenated as above
  renorm = function(x) x |> NormalizeData() |> ScaleData(features = rownames(x))
  human_p = renorm(human[!rownames(human) %in% drop, ])
  mouse_p = renorm(mouse_subset[!rownames(mouse_subset) %in% hs_to_mm[drop], ])
  stopifnot(nrow(human) - nrow(human_p) == length(drop),
            nrow(mouse_subset) - nrow(mouse_p) == length(drop))
  concat_p = prepare_concat(concatenate(human_p, mouse_p, homology_tbl))
  stopifnot(setequal(rownames(concat_p), setdiff(rownames(concat), drop)))

  out_dir = sprintf('data/genes%d/', round(10 * frac))
  save_concat(concat_p, paste0(out_dir, '03_concatenated_data'))
  save_species_inputs(concat_p, human_p, mouse_p, out_dir)
  print(sprintf('Dropped %d of %d shared genes -> %s (cells kept: %d of %d)',
                length(drop), n_shared, out_dir, ncol(concat_p), ncol(concat)))
}
