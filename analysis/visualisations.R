# Figures for the cross-species integration benchmark. Run from the repository root.
# Each tissue is described once in `tissues` (below); every figure is built from that.
# Output layout:
#   figures/<tissue>/scores/scores_<case>.png     batch vs bio scatter, one per case
#   figures/<tissue>/scores/pareto_<case>.png     same, with the Pareto front labelled
#   figures/<tissue>/scores/scores_all_cases.png  all cases, shared axes
#   figures/<tissue>/scores/scores.csv            the plotted scores
#   figures/<tissue>/umap/<case>/<NN>_<method>.png  species | cell type, evaluation cells
#   figures/<tissue>/spatial/<sample>/<annotator>.png  ImageDimPlot coloured by each annotator
#   figures/<tissue>/annotation/marker_auroc.csv, marker_auroc_<species>.png  how well each annotator's labels match panel markers
source('library.R')
source('source.R')
library(rhdf5)

# ---- tissues and their cases ----
# Each case: concatenated (unintegrated) data, 04 integration outputs, 05 scores.
# Add the other tissues here when their pipelines use the same layout.
bone_marrow_cases = list(
  full      = list(concat = 'bone_marrow/data/03_concatenated_data/03_concatenated_data_full',
                   int    = 'bone_marrow/data/04_integrated_data/full_mouse/',
                   scores = 'bone_marrow/data/05_evaluated_score/full_mouse/'),
  dissected = list(concat = 'bone_marrow/data/03_concatenated_data/03_concatenated_data_subset',
                   int    = 'bone_marrow/data/04_integrated_data/dissected_mouse/',
                   scores = 'bone_marrow/data/05_evaluated_score/dissected_mouse/')
)
for (g in c(2, 4, 6, 8)) {
  bone_marrow_cases[[paste0('genes', g)]] = list(
    concat = str_glue('bone_marrow/data/genes{g}/03_concatenated_data'),
    int    = str_glue('bone_marrow/data/genes{g}/'),
    scores = str_glue('bone_marrow/data/genes{g}/'))
}
# annotated single-species objects (02) for the spatial plots
bone_marrow_spatial = list(
  human           = 'bone_marrow/data/02_annotated_data/02_human_annotated.rds',
  mouse_full      = 'bone_marrow/data/02_annotated_data/02_mouse_annotated.rds',
  mouse_dissected = 'bone_marrow/data/02_annotated_data/02_mouse_subset_annotated.rds'
)
# annotation quality: annotated object + the 02 script whose `markers = list(...)` to use
bone_marrow_annotation = list(
  human = list(obj = 'bone_marrow/data/02_annotated_data/02_human_annotated.rds',
               markers = 'bone_marrow/02_human_annotate.R'),
  mouse = list(obj = 'bone_marrow/data/02_annotated_data/02_mouse_annotated.rds',
               markers = 'bone_marrow/02_mouse_annotate.R')
)
colon_cases = list(
  full      = list(concat = 'colon/data/03_concatenated_data/03_concatenated_data_full',
                   int    = 'colon/data/04_integrated_data/full_mouse/',
                   scores = 'colon/data/05_evaluated_score/full_mouse/'),
  dissected = list(concat = 'colon/data/03_concatenated_data/03_concatenated_data_subset',
                   int    = 'colon/data/04_integrated_data/dissected_mouse/',
                   scores = 'colon/data/05_evaluated_score/dissected_mouse/')
)
for (g in c(2, 4, 6, 8)) {
  colon_cases[[paste0('genes', g)]] = list(
    concat = str_glue('colon/data/genes{g}/03_concatenated_data'),
    int    = str_glue('colon/data/genes{g}/'),
    scores = str_glue('colon/data/genes{g}/'))
}
colon_spatial = list(
  human           = 'colon/data/02_annotated_data/02_human_annotated.rds',
  mouse_full      = 'colon/data/02_annotated_data/02_mouse_annotated.rds',
  mouse_dissected = 'colon/data/02_annotated_data/02_mouse_subset_annotated.rds'
)
colon_annotation = list(
  human = list(obj = 'colon/data/02_annotated_data/02_human_annotated.rds',
               markers = 'colon/02_human_annotated.R'),
  mouse = list(obj = 'colon/data/02_annotated_data/02_mouse_annotated.rds',
               markers = 'colon/02_mouse_annotate.qmd')
)
fetal_liver_cases = list(
  full      = list(concat = 'fetal_liver/data/03_concatenated_data/03_concatenated_data_full',
                   int    = 'fetal_liver/data/04_integrated_data/full_mouse/',
                   scores = 'fetal_liver/data/05_evaluated_score/full_mouse/'),
  dissected = list(concat = 'fetal_liver/data/03_concatenated_data/03_concatenated_data_subset',
                   int    = 'fetal_liver/data/04_integrated_data/dissected_mouse/',
                   scores = 'fetal_liver/data/05_evaluated_score/dissected_mouse/')
)
for (g in c(2, 4, 6, 8)) {
  fetal_liver_cases[[paste0('genes', g)]] = list(
    concat = str_glue('fetal_liver/data/genes{g}/03_concatenated_data'),
    int    = str_glue('fetal_liver/data/genes{g}/'),
    scores = str_glue('fetal_liver/data/genes{g}/'))
}
fetal_liver_spatial = list(
  human           = 'fetal_liver/data/02_annotated_data/02_human_annotated.rds',
  mouse_full      = 'fetal_liver/data/02_annotated_data/02_mouse_annotated.rds',
  mouse_dissected = 'fetal_liver/data/02_annotated_data/02_mouse_subset_annotated.rds'
)
fetal_liver_annotation = list(
  human = list(obj = 'fetal_liver/data/02_annotated_data/02_human_annotated.rds',
               markers = 'fetal_liver/02_human_annotate.R'),
  mouse = list(obj = 'fetal_liver/data/02_annotated_data/02_mouse_annotated.rds',
               markers = 'fetal_liver/02_mouse_annotate.R')
)
breast_cancer_cases = list()
for (d in c('full_mouse', 'dissected_mouse', 'genes2', 'genes4', 'genes6', 'genes8')) {
  breast_cancer_cases[[sub('_mouse$', '', d)]] = list(
    concat = str_glue('breast_cancer/data/{d}/03_concatenated_data'),
    int    = str_glue('breast_cancer/data/{d}/'),
    scores = str_glue('breast_cancer/data/{d}/'))
}
breast_cancer_spatial = list(
  human           = 'breast_cancer/data/02_annotated_data/02_human_annotated.rds',
  mouse_full      = 'breast_cancer/data/02_annotated_data/02_mouse_annotated.rds',
  mouse_dissected = 'breast_cancer/data/dissected_mouse/02_mouse_annotated.rds'
)
# breast_cancer's marker-group names differ from its labels: label_map gives label -> marker group;
# the mouse notebook's list is called markers_pymt
breast_cancer_annotation = list(
  human = list(obj = 'breast_cancer/data/02_annotated_data/02_human_annotated.rds',
               markers = 'breast_cancer/02_human_annotate.qmd',
               label_map = c('B-cells' = 'B', 'CAFs' = 'CAF', 'Cancer Epithelial' = 'Tumour',
                             'Endothelial' = 'Endothelial', 'Myeloid' = 'Myeloid',
                             'Normal Epithelial' = 'Basal', 'Plasmablasts' = 'Plasmablast',
                             'PVL' = 'PVL', 'T-cells' = 'T_NK')),
  mouse = list(obj = 'breast_cancer/data/02_annotated_data/02_mouse_annotated.rds',
               markers = 'breast_cancer/02_mouse_annotate.qmd', markers_var = 'markers_pymt',
               label_map = c('B-cells' = 'B_cells', 'CAFs' = 'CAF', 'Cancer Epithelial' = 'Tumour',
                             'Endothelial' = 'Endothelial', 'Myeloid' = 'Myeloid',
                             'Normal Epithelial' = 'Basal', 'Plasmablasts' = 'Plasma',
                             'PVL' = 'PVL', 'T-cells' = 'T_cells'))
)
tissues = list(bone_marrow = list(cases = bone_marrow_cases, spatial = bone_marrow_spatial,
                                  annotation = bone_marrow_annotation),
               colon       = list(cases = colon_cases, spatial = colon_spatial,
                                  annotation = colon_annotation),
               fetal_liver = list(cases = fetal_liver_cases, spatial = fetal_liver_spatial,
                                  annotation = fetal_liver_annotation),
               breast_cancer = list(cases = breast_cancer_cases, spatial = breast_cancer_spatial,
                                    annotation = breast_cancer_annotation))
# tissue names on the command line run only those, e.g. Rscript analysis/visualisations.R colon
run_only = commandArgs(TRUE)
if (length(run_only)) {
  stopifnot(all(run_only %in% names(tissues)))
  tissues = tissues[run_only]
}

# ---- methods ----
# h5ad file in the case's 04 folder and the obsm key holding its own UMAP
# (a plain X_umap in most 04 files is the unintegrated UMAP carried over)
methods = tribble(
  ~method,              ~file,                        ~umap_key,
  'Unintegrated',       NA,                           'X_umap',
  'Unintegrated UINMF', NA,                           NA,          # computed from the UINMF H matrix
  'FastMNN',            '04.01_FastMNN.h5ad',         'X_umap_mnn',
  'SAMap',              '04.02_SAMap.h5ad',           'X_umap',
  'CCA',                '04.03_CCA.h5ad',             'X_umap_cca',
  'Harmony',            '04.04_Harmony.h5ad',         'X_umap_harmony',
  'Scanorama',          '04.05_Scanorama.h5ad',       'X_umap_scanorama',
  'RPCA',               '04.06_RPCA.h5ad',            'X_umap_rpca',
  'scVI',               '04.07_scVI.h5ad',            'X_umap_scvi',
  'scANVI',             '04.08_scANVI.h5ad',          'X_umap_scanvi',
  'BBKNN',              '04.09_BBKNN.h5ad',           'X_umap_bbknn',
  'LIGER',              '04.10_LIGER.h5ad',           'X_umap_liger',
  'LIGER UINMF',        '04.11_LIGER_UINMF.h5ad',     'X_umap',
  'FastMNN UINMF',      '04.12_FastMNN_UINMF.h5ad',   'X_umap_mnn_uinmf',
  'CCA UINMF',          '04.13_CCA_UINMF.h5ad',       'X_umap_cca_uinmf',
  'Harmony UINMF',      '04.14_Harmony_UINMF.h5ad',   'X_umap_harmony_uinmf',
  'Scanorama UINMF',    '04.15_Scanorama_UINMF.h5ad', 'X_umap_scanorama_uinmf',
  'RPCA UINMF',         '04.16_RPCA_UINMF.h5ad',      'X_umap_rpca_uinmf',
  'BBKNN UINMF',        '04.17_BBKNN_UINMF.h5ad',     'X_umap_bbknn_uinmf',
  # the same methods on an scVI latent trained without the species key (04.18), instead of PCA
  'Unintegrated scVI-latent', '04.18_scVI_latent.h5ad', 'X_umap_scvi_latent',
  'FastMNN scVI-latent',   '04.19_FastMNN_scVI.h5ad',   'X_umap_mnn_scvi',
  'CCA scVI-latent',       '04.20_CCA_scVI.h5ad',       'X_umap_cca_scvi',
  'Harmony scVI-latent',   '04.21_Harmony_scVI.h5ad',   'X_umap_harmony_scvi',
  'Scanorama scVI-latent', '04.22_Scanorama_scVI.h5ad', 'X_umap_scanorama_scvi',
  'RPCA scVI-latent',      '04.23_RPCA_scVI.h5ad',      'X_umap_rpca_scvi',
  'BBKNN scVI-latent',     '04.24_BBKNN_scVI.h5ad',     'X_umap_bbknn_scvi'
) |>
  mutate(base  = str_remove(method, ' (UINMF|scVI-latent)$'),
         # input the method runs on: PCA of shared genes, UINMF H matrix, or scVI latent
         uinmf = case_when(str_detect(method, 'UINMF')       ~ 'UINMF',
                           str_detect(method, 'scVI-latent') ~ 'scVI-latent',
                           TRUE                              ~ 'Base'),
         # file-name prefix: 04.xx number, 00/00b for the two unintegrated versions
         id    = coalesce(str_extract(file, '^04\\.\\d+'),
                          if_else(method == 'Unintegrated', '04.00', '04.00b')))

base_cols = c(Unintegrated = '#000000',
              FastMNN      = '#E69F00',
              SAMap        = '#0072B2',
              CCA          = '#009E73',
              Scanorama    = '#D55E00',
              RPCA         = '#332288',
              Harmony      = '#AA4499',
              BBKNN        = '#999933',
              scVI         = '#882255',
              scANVI       = '#CC6677',
              LIGER        = '#8C510A')
species_cols = c(hsapiens = '#4477AA', mmusculus = '#EE6677')

save_plot = function(p, file, width = 11, height = 7) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggsave(file, p, width = width, height = height, dpi = 300, bg = 'white')
}

# ---- scores: batch removal vs biological conservation ----
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label',
                  'isolated_label_F1', 'isolated_label_silhouette', 'cLISI')

# one 05 csv (rows = metrics, one column named after the method) -> long format
read_metrics = function(f) {
  df = read_csv(f, show_col_types = FALSE)
  tibble(metric = df[[1]], value = as.numeric(df[[2]]), method = names(df)[2])
}

# per method and case: mean of the available batch / bio metrics
case_scores = function(cases) {
  imap_dfr(cases, function(cs, case) {
    files = list.files(cs$scores, pattern = '^05\\..*\\.csv$', full.names = TRUE)
    map_dfr(files, read_metrics) |>
      mutate(group = case_when(metric %in% batch_metrics ~ 'batch',
                               metric %in% bio_metrics   ~ 'bio')) |>
      filter(!is.na(group), !is.na(value)) |>
      group_by(method, group) |>
      summarise(score = mean(value), .groups = 'drop') |>
      pivot_wider(names_from = group, values_from = score) |>
      mutate(case = case)
  }) |>
    left_join(select(methods, method, base, uinmf), by = 'method') |>
    mutate(base  = factor(base, levels = names(base_cols)),
           uinmf = factor(uinmf, levels = c('Base', 'UINMF', 'scVI-latent')),
           case  = factor(case, levels = names(cases)))
}

# colour = base algorithm, shape = gene set, dashed lines = Unintegrated's scores
plot_scores = function(s, title = NULL) {
  ref = s |> filter(method == 'Unintegrated')
  ggplot(s, aes(batch, bio)) +
    geom_vline(data = ref, aes(xintercept = batch, linetype = 'Unintegrated'), colour = 'grey40') +
    geom_hline(data = ref, aes(yintercept = bio,   linetype = 'Unintegrated'), colour = 'grey40') +
    geom_point(aes(colour = base, shape = uinmf), size = 3.5, stroke = 1.2) +
    scale_colour_manual(values = base_cols, name = 'Base algorithm', drop = FALSE) +
    scale_shape_manual(values = c(Base = 16, UINMF = 2, `scVI-latent` = 0), name = 'Input', drop = FALSE,
                       labels = c(Base = 'Shared orthologs (PCA)', UINMF = 'Shared + unshared (UINMF)',
                                  `scVI-latent` = 'scVI latent (shared orthologs)')) +
    scale_linetype_manual(values = c(Unintegrated = 'dashed'), name = NULL) +
    guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2),
           linetype = guide_legend(order = 3)) +
    labs(title = title, x = 'Batch removal score', y = 'Biological conservation score') +
    theme_classic(base_size = 14) +
    theme(strip.background = element_blank())
}

# Pareto front: a method is on it if no other method in the same case is at least
# as good on both scores and strictly better on one. Front points are drawn full
# size and labelled, the rest dimmed; the line is the staircase edge of the
# region the front dominates (across to the next front point, then down).
pareto_flag = function(batch, bio) {
  map2_lgl(batch, bio, \(b, s) !any(batch >= b & bio >= s & (batch > b | bio > s)))
}

plot_pareto = function(s, title = NULL) {
  s = s |> mutate(pareto = pareto_flag(batch, bio))
  front = s |> filter(pareto) |> arrange(batch)
  ref = s |> filter(method == 'Unintegrated')
  ggplot(s, aes(batch, bio)) +
    geom_vline(data = ref, aes(xintercept = batch, linetype = 'Unintegrated'), colour = 'grey40') +
    geom_hline(data = ref, aes(yintercept = bio,   linetype = 'Unintegrated'), colour = 'grey40') +
    geom_step(data = front, direction = 'hv', colour = 'grey30', linewidth = 0.7) +
    geom_point(data = filter(s, !pareto), aes(colour = base, shape = uinmf),
               size = 3, stroke = 1, alpha = 0.25) +
    geom_point(data = front, aes(colour = base, shape = uinmf), size = 4.5, stroke = 1.6) +
    ggrepel::geom_text_repel(data = front, aes(label = method, colour = base),
                             size = 4.5, fontface = 'bold', show.legend = FALSE,
                             max.overlaps = Inf, min.segment.length = 0,
                             box.padding = 0.6, seed = 1) +
    scale_colour_manual(values = base_cols, name = 'Base algorithm', drop = FALSE) +
    scale_shape_manual(values = c(Base = 16, UINMF = 2, `scVI-latent` = 0), name = 'Input', drop = FALSE,
                       labels = c(Base = 'Shared orthologs (PCA)', UINMF = 'Shared + unshared (UINMF)',
                                  `scVI-latent` = 'scVI latent (shared orthologs)')) +
    scale_linetype_manual(values = c(Unintegrated = 'dashed'), name = NULL) +
    guides(colour = guide_legend(order = 1, override.aes = list(alpha = 1, size = 3.5, stroke = 1.2)),
           shape = guide_legend(order = 2, override.aes = list(alpha = 1, size = 3.5, stroke = 1.2)),
           linetype = guide_legend(order = 3)) +
    labs(title = title, x = 'Batch removal score', y = 'Biological conservation score',
         caption = 'Line: Pareto front. Labelled points are Pareto-optimal; dimmed points are dominated.') +
    theme_classic(base_size = 14) +
    theme(plot.caption = element_text(hjust = 0, size = 10))
}

# ---- UMAPs ----
# obs names and one obsm embedding from an h5ad, without loading the whole file
read_h5ad_embedding = function(file, key) {
  idx = h5readAttributes(file, 'obs')[['_index']]
  cells = as.character(h5read(file, paste0('obs/', idx)))
  emb = t(h5read(file, paste0('obsm/', key)))[, 1:2]
  h5closeAll()
  data.frame(cell = cells, umap_1 = emb[, 1], umap_2 = emb[, 2])
}

# match a method's cell names to the concatenated data's: SAMap drops the
# species prefix (barcodes do not overlap between species, so the bare barcode is unique)
match_cells = function(cells, reference) {
  i = match(cells, reference)
  bare = str_remove(reference, '^(hsapiens|mmusculus)_')
  i[is.na(i)] = match(cells[is.na(i)], bare)
  reference[i]
}

# UMAP of the UINMF H matrix (no integration on top), as RunUMAP would by default
uinmf_h_umap = function(int_dir, cells) {
  H = read.csv(file.path(int_dir, 'uinmf_H_raw.csv'), row.names = 1, check.names = FALSE)
  H = as.matrix(H[intersect(cells, rownames(H)), ])
  emb = uwot::umap(H, n_neighbors = 30, metric = 'cosine', min_dist = 0.3, seed = 42)
  data.frame(cell = rownames(H), umap_1 = emb[, 1], umap_2 = emb[, 2])
}

plot_umap = function(emb, meta, title, celltype_cols) {
  d = emb |>
    inner_join(meta, by = 'cell') |>
    slice_sample(prop = 1)                     # random draw order: no species on top
  # zoom to the central 99% of cells on each axis; a few scattered outliers
  # would otherwise shrink the main structure (the subtitle reports how many)
  xl = quantile(d$umap_1, c(0.005, 0.995)); yl = quantile(d$umap_2, c(0.005, 0.995))
  xl = xl + c(-0.05, 0.05) * diff(xl);      yl = yl + c(-0.05, 0.05) * diff(yl)
  n_out = sum(d$umap_1 < xl[1] | d$umap_1 > xl[2] | d$umap_2 < yl[1] | d$umap_2 > yl[2])
  panel = function(col, values, name) {
    ggplot(d, aes(umap_1, umap_2, colour = .data[[col]])) +
      geom_point(size = 0.3, stroke = 0) +
      scale_colour_manual(values = values, name = name) +
      guides(colour = guide_legend(override.aes = list(size = 3))) +
      coord_equal(xlim = xl, ylim = yl) +
      theme_void(base_size = 12)
  }
  (panel('species', species_cols, 'Species') | panel('celltype', celltype_cols, 'Cell type')) +
    patchwork::plot_annotation(title = title,
                               subtitle = sprintf('%d cells (%d outside the view)', nrow(d), n_out))
}

# ---- spatial (ImageDimPlot) by annotator ----
annotators = c(seurat        = 'celltype_seurat',          # same-species reference
               singleR       = 'celltype_singleR',
               seurat_cross  = 'celltype_seurat_cross',    # other species' reference
               singleR_cross = 'celltype_singleR_cross',
               ebcc          = 'ebcc_label')

spatial_sample = function(file, sample, out_dir, overwrite = FALSE) {
  outs = file.path(out_dir, sample, paste0(names(annotators), '.png'))
  if (all(file.exists(outs)) && !overwrite) return(invisible())
  obj = readRDS(file)
  # one colour per label, shared by every annotator; SingleR's pruned (NA) calls are grey
  for (a in annotators) obj@meta.data[[a]] = replace_na(as.character(obj@meta.data[[a]]), 'Unassigned')
  labs = sort(setdiff(unique(unlist(obj@meta.data[annotators])), 'Unassigned'))
  cols = c(setNames(scales::hue_pal()(length(labs)), labs), Unassigned = 'grey50')
  for (k in seq_along(annotators)) {
    # (border.color = NA hides the centroid points entirely, so keep the default border)
    p = ImageDimPlot(obj, fov = Images(obj)[1], group.by = annotators[[k]], cols = cols,
                     size = 0.8, axes = TRUE) +
      labs(title = sprintf('%s - %s', sample, names(annotators)[k]))
    save_plot(p, outs[k], width = 10, height = 9)
    message('saved ', outs[k])
  }
}

# ---- annotation quality: marker AUROC per annotator ----
# For each cell type: score every cell by the mean log-normalised expression of
# that type's panel markers; AUROC = P(a cell labelled with the type scores
# higher than a cell labelled otherwise). 0.5 = label unrelated to its markers.
read_markers = function(script, var = 'markers') {
  l = readLines(script)
  s = grep(paste0('^', var, ' = list'), l)[1]
  e = s - 1 + which(l[s:length(l)] == ')')[1]
  env = new.env()
  eval(parse(text = l[s:e]), envir = env)
  env[[var]]
}

# a tissue's annotation entry -> markers named by label (optional markers_var, label_map = label -> group)
annotation_markers = function(a) {
  m = read_markers(a$markers, if (is.null(a$markers_var)) 'markers' else a$markers_var)
  if (!is.null(a$label_map)) m = setNames(m[a$label_map], names(a$label_map))
  m[!vapply(m, is.null, logical(1))]
}

auroc = function(score, pos) {            # Mann-Whitney form, ties averaged
  r = rank(score); n1 = as.numeric(sum(pos)); n0 = as.numeric(sum(!pos))   # numeric: n1 * n0 overflows integer
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

marker_auroc = function(obj_file, markers) {
  obj = readRDS(obj_file)
  assay = DefaultAssay(obj)
  rownames(obj) = obj[[assay]]@meta.data$gene_name
  X = LayerData(obj, assay = assay, layer = 'data')
  markers = lapply(markers, intersect, rownames(X))
  meta = obj@meta.data
  rm(obj); gc()
  # labels and marker-group names are matched ignoring case and punctuation
  # (e.g. colon mouse markers 'Smooth_Muscle' for the label 'Smooth muscle')
  key = function(x) tolower(gsub('[^A-Za-z0-9]+', '_', x))
  imap_dfr(annotators, function(col, annotator) {
    lab = key(as.character(meta[[col]]))
    ok = !is.na(lab)
    called = unique(lab[ok])
    map_dfr(names(markers), function(k) {
      tibble(annotator = annotator, celltype = k, n_markers = length(markers[[k]]),
             n_cells = sum(lab[ok] == key(k)),
             auroc = if (length(markers[[k]]) == 0 || !key(k) %in% called) NA_real_ else
               auroc(colMeans(as.matrix(X[markers[[k]], ok, drop = FALSE])), lab[ok] == key(k)))
    })
  })
}

# one species' marker AUROC heatmap; the summary row is each annotator's mean over
# the cell types it calls (annotators may cover different types)
plot_marker_auroc = function(res, title = NULL) {
  # drop marker groups no annotator calls (e.g. mouse Macrophage/Endothelial: markers only)
  res = res |> group_by(celltype) |> filter(any(!is.na(auroc))) |> ungroup()
  mean_row = 'Mean (all called types)'
  means = res |> group_by(annotator) |>
    summarise(auroc = mean(auroc, na.rm = TRUE)) |> mutate(celltype = mean_row)
  is_mean = function(x) x == mean_row
  d = bind_rows(res, means) |>
    mutate(annotator = factor(annotator, levels = names(annotators)),
           row = if_else(is_mean(celltype), celltype, sprintf('%s (%d)', celltype, n_markers)))
  # rows: cell types in marker order, the mean at the bottom
  d$row = factor(d$row, levels = rev(unique(c(d$row[!is_mean(d$celltype)], mean_row))))
  ggplot(d, aes(annotator, row, fill = auroc)) +
    geom_tile(colour = 'white') +
    geom_hline(yintercept = 1.5, linewidth = 0.8) +
    geom_text(aes(label = ifelse(is.na(auroc), '–', sprintf('%.2f', auroc)),
                  fontface = ifelse(is_mean(celltype), 'bold', 'plain')), size = 4) +
    scale_fill_gradient2(low = '#B2182B', mid = 'white', high = '#2166AC', midpoint = 0.5,
                         limits = c(0.3, 1), oob = scales::squish, na.value = 'grey90',
                         name = 'Marker AUROC') +
    labs(title = title, x = NULL, y = 'Cell type (number of panel markers)',
         caption = paste('AUROC of the mean marker expression for cells given each label vs all other cells.',
                         '0.5 = no better than chance; "–" = no panel markers or label never called.',
                         'seurat/singleR use the same-species reference, *_cross the other species\'.',
                         'Mean: over the cell types each annotator calls.',
                         sep = '\n')) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(),
          plot.caption = element_text(hjust = 0, size = 9), plot.title = element_text(face = 'bold'))
}

# ---- annotation agreement: seurat vs singleR (same-species reference) ----
# per species: share of the cells Seurat gives each label that SingleR gives each
# label (rows sum to 1); SingleR's pruned (NA) calls are 'Unassigned'
annotation_agreement = function(obj_file) {
  meta = readRDS(obj_file)@meta.data; gc()
  tibble(seurat  = replace_na(as.character(meta$celltype_seurat), 'Unassigned'),
         singleR = replace_na(as.character(meta$celltype_singleR), 'Unassigned')) |>
    dplyr::count(seurat, singleR) |>
    group_by(seurat) |> mutate(frac = n / sum(n)) |> ungroup()
}

# base_size / text_size shrink it for the combined cross-tissue figure
plot_annotation_agreement = function(res, title = NULL, base_size = 12, text_size = 3.3,
                                     axis_titles = c('SingleR label (% of cells)', 'Seurat label (% of cells)')) {
  # same label order on both axes so agreement sits on the diagonal; Unassigned last
  labs = sort(setdiff(union(res$seurat, res$singleR), 'Unassigned'))
  labs = c(labs, intersect('Unassigned', c(res$seurat, res$singleR)))
  row_n = res |> group_by(seurat) |> summarise(n = sum(n))
  row_pct = 100 * replace_na(row_n$n[match(labs, row_n$seurat)], 0L) / sum(res$n)   # share of all cells given that Seurat label
  pct_fmt = \(x) ifelse(x > 0 & x < 0.05, '<0.1%', sprintf('%.1f%%', x))
  row_lab = setNames(sprintf('%s (%s)', labs, pct_fmt(row_pct)), labs)
  col_n = res |> group_by(singleR) |> summarise(n = sum(n))
  col_pct = 100 * replace_na(col_n$n[match(labs, col_n$singleR)], 0L) / sum(res$n)   # share of all cells given that SingleR label
  col_lab = setNames(sprintf('%s (%s)', labs, pct_fmt(col_pct)), labs)
  agree = sum(res$n[res$seurat == res$singleR]) / sum(res$n)
  # Cohen's kappa: agreement beyond what the two label frequencies give by chance (Unassigned is a category)
  p_seurat = res |> group_by(lab = seurat) |> summarise(p = sum(n)) |> mutate(p = p / sum(p))
  p_singleR = res |> group_by(lab = singleR) |> summarise(p = sum(n)) |> mutate(p = p / sum(p))
  chance = sum(inner_join(p_seurat, p_singleR, by = 'lab') |> with(p.x * p.y))
  kappa = (agree - chance) / (1 - chance)
  d = res |> select(-species) |>
    complete(seurat = labs, singleR = labs, fill = list(n = 0)) |>
    # labels Seurat never calls have no row fractions (grey)
    group_by(seurat) |> mutate(frac = if (all(n == 0)) NA else coalesce(frac, 0)) |> ungroup() |>
    filter(seurat != 'Unassigned') |>         # only SingleR prunes, so this row is always empty
    mutate(diag    = seurat == singleR,
           seurat  = factor(row_lab[seurat], levels = rev(row_lab[setdiff(labs, 'Unassigned')])),
           singleR = factor(singleR, levels = labs))
  ggplot(d, aes(singleR, seurat, fill = frac)) +
    geom_tile(colour = 'white') +
    geom_tile(data = filter(d, diag), colour = 'black', fill = NA, linewidth = 0.6) +
    geom_text(aes(label = ifelse(is.na(frac) | frac < 0.005, '', sprintf('%.0f', 100 * frac)),
                  colour = !is.na(frac) & frac > 0.6), size = text_size, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = 'grey15', `TRUE` = 'white')) +
    scale_fill_gradient(low = '#F7FBFF', high = '#08306B', limits = c(0, 1), na.value = 'grey90',
                        labels = scales::percent, name = '% of Seurat\nlabel') +
    scale_x_discrete(labels = col_lab) +
    labs(title = title, subtitle = sprintf('Overall agreement: %.1f%% of %s cells; Cohen\'s \u03ba = %.2f',
                                           100 * agree, scales::comma(sum(res$n)), kappa),
         x = axis_titles[1], y = axis_titles[2]) +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(),
          plot.title = element_text(face = 'bold'))
}

# ---- annotation quality on the evaluation half: annotators, majority vote, EBCC ----
# Computed by analysis/vote_vs_ebcc.R (all labellings scored on the same cells, with bootstrap
# intervals), which also draws heldout_auroc*.png; the label sets, the vote and the
# per-cell-type heatmap are defined here.
label_sets = c('Seurat', 'SingleR', 'Seurat (cross)', 'SingleR (cross)',
               'Majority vote', 'EBCC (f = 0)', 'EBCC (chosen f)')

# most frequent label over the four annotators (NA calls skipped); ties -> NA
majority_vote = function(meta) {
  L = sapply(annotators[c('seurat', 'singleR', 'seurat_cross', 'singleR_cross')],
             \(col) as.character(meta[[col]]))
  apply(L, 1, function(x) {
    tab = table(x[!is.na(x)])
    top = names(tab)[tab == max(tab)]
    if (length(top) == 1) top else NA_character_
  })
}

# per-class rows -> each label set's mean over the types it calls
heldout_means = function(res) {
  res |> group_by(across(any_of(c('tissue', 'species', 'chosen_frac'))), label_set) |>
    summarise(mean_auroc = mean(auroc, na.rm = TRUE), n_types = sum(!is.na(auroc)), .groups = 'drop')
}

# one tissue: rows = cell types, columns = label sets, one facet per species; grey
# rows = no markers on that species' panel; '-' = label set never calls the type
plot_heldout_celltype = function(res, title = NULL, base_size = 10, text_size = 2.6) {
  key = function(x) tolower(gsub('[^A-Za-z0-9]+', '_', x))
  res = filter(res, label_set %in% label_sets) |>
    # one name per type across species (colon mouse markers spell 'Smooth_Muscle'): first seen
    mutate(k = key(celltype), celltype = celltype[match(k, k)]) |>
    # drop types no labelling calls in either species (e.g. bone_marrow Macrophage: markers only)
    group_by(k) |> filter(any(n_cells > 0)) |> ungroup() |> select(-k)
  mean_row = 'Mean (called types)'
  d = bind_rows(res, heldout_means(res) |> rename(auroc = mean_auroc) |>
                  mutate(celltype = mean_row, n_markers = NA)) |>
    # types that do not exist in one species get a blank tile
    complete(species, celltype, label_set) |>
    group_by(species, celltype) |>
    mutate(absent = all(is.na(n_cells)) & celltype != mean_row,
           no_markers = !absent & all(n_markers %in% 0)) |> ungroup() |>
    mutate(label_set = factor(label_set, levels = label_sets),
           celltype  = factor(celltype, levels = rev(c(unique(res$celltype), mean_row))),
           text = case_when(no_markers | absent ~ '', is.na(auroc) ~ '–', TRUE ~ sprintf('%.2f', auroc)))
  ggplot(d, aes(label_set, celltype, fill = auroc)) +
    geom_tile(colour = 'white') +
    geom_tile(data = filter(d, no_markers), fill = 'grey65', colour = 'white') +
    geom_tile(data = filter(d, absent), fill = 'white', colour = 'white') +
    geom_text(aes(label = text, fontface = ifelse(celltype == mean_row, 'bold', 'plain')), size = text_size) +
    geom_hline(yintercept = 1.5, linewidth = 0.6) +
    geom_vline(xintercept = 5.5, linewidth = 0.6) +                     # EBCC to the right
    facet_wrap(~ species, nrow = 1) +
    scale_fill_gradient2(low = '#B2182B', mid = 'white', high = '#2166AC', midpoint = 0.5,
                         limits = c(0.3, 1), oob = scales::squish, na.value = 'grey92',
                         name = 'Marker AUROC') +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(),
          strip.text = element_text(face = 'bold', hjust = 0), plot.title = element_text(face = 'bold'))
}

umap_case = function(cs, case, out_dir, overwrite = FALSE) {
  concat = readRDS(paste0(cs$concat, '.rds'))
  meta = data.frame(cell     = colnames(concat),
                    species  = concat$orig.ident,
                    celltype = as.character(concat$celltype))
  rm(concat); gc()
  eval_cells = readLines(file.path(cs$int, 'eval_cells.txt'))
  meta = meta[meta$cell %in% eval_cells, ]
  celltype_cols = setNames(scales::hue_pal()(length(unique(meta$celltype))),
                           sort(unique(meta$celltype)))

  for (k in seq_len(nrow(methods))) {
    m = methods[k, ]
    out = file.path(out_dir, case, sprintf('%s_%s.png', m$id, str_replace_all(m$method, ' ', '_')))
    if (file.exists(out) && !overwrite) next
    # methods this tissue has not run (e.g. the scVI-latent variants) are skipped
    if (!is.na(m$file) && !file.exists(file.path(cs$int, m$file))) next
    emb = if (m$method == 'Unintegrated') {
      read_h5ad_embedding(paste0(cs$concat, '.h5ad'), m$umap_key)
    } else if (m$method == 'Unintegrated UINMF') {
      uinmf_h_umap(cs$int, meta$cell)
    } else {
      read_h5ad_embedding(file.path(cs$int, m$file), m$umap_key)
    }
    emb$cell = match_cells(emb$cell, meta$cell)
    p = plot_umap(emb[!is.na(emb$cell), ], meta, sprintf('%s - %s', m$method, case), celltype_cols)
    save_plot(p, out, width = 14, height = 6)
    message('saved ', out)
  }
}

# ---- run ----
for (tissue in names(tissues)) {
  cases = tissues[[tissue]]$cases
  fig_dir = file.path('figures', tissue)

  for (sample in names(tissues[[tissue]]$spatial)) {
    spatial_sample(tissues[[tissue]]$spatial[[sample]], sample, file.path(fig_dir, 'spatial'))
  }

  ann = tissues[[tissue]]$annotation
  auroc_res = imap_dfr(ann, \(a, sp) marker_auroc(a$obj, annotation_markers(a)) |> mutate(species = sp))
  dir.create(file.path(fig_dir, 'annotation'), showWarnings = FALSE, recursive = TRUE)
  write_csv(auroc_res, file.path(fig_dir, 'annotation', 'marker_auroc.csv'))
  for (sp in names(ann)) {
    save_plot(plot_marker_auroc(filter(auroc_res, species == sp), title = str_glue('{tissue} - {sp}')),
              file.path(fig_dir, 'annotation', str_glue('marker_auroc_{sp}.png')), width = 8.5, height = 8)
  }
  agree_res = imap_dfr(ann, \(a, sp) annotation_agreement(a$obj) |> mutate(species = sp))
  write_csv(agree_res, file.path(fig_dir, 'annotation', 'seurat_singleR_agreement.csv'))
  p = purrr::map(names(ann), \(sp) plot_annotation_agreement(filter(agree_res, species == sp), title = sp))
  save_plot(patchwork::wrap_plots(p, nrow = 1) +
              patchwork::plot_annotation(
                title = str_glue('{tissue}: Seurat vs SingleR, same-species reference'),
                caption = paste('Cell values: % of cells with that Seurat label given each SingleR label (rows sum to 100).',
                                'Outlined: both annotators agree. Unassigned = SingleR pruned call.'),
                theme = theme(plot.title = element_text(face = 'bold', size = 15))),
            file.path(fig_dir, 'annotation', 'seurat_singleR_agreement.png'), width = 18, height = 8.5)

  scores = case_scores(cases)
  dir.create(file.path(fig_dir, 'scores'), showWarnings = FALSE, recursive = TRUE)
  write_csv(scores, file.path(fig_dir, 'scores', 'scores.csv'))
  for (case in names(cases)) {
    s = filter(scores, case == !!case)
    save_plot(plot_scores(s, title = str_glue('{tissue} - {case}')),
              file.path(fig_dir, 'scores', str_glue('scores_{case}.png')))
    save_plot(plot_pareto(s, title = str_glue('{tissue} - {case}')),
              file.path(fig_dir, 'scores', str_glue('pareto_{case}.png')))
  }
  save_plot(plot_scores(scores, title = tissue) + facet_wrap(~ case, nrow = 2),
            file.path(fig_dir, 'scores', 'scores_all_cases.png'), width = 15, height = 9)

  for (case in names(cases)) umap_case(cases[[case]], case, file.path(fig_dir, 'umap'))
}

# ---- all tissues: seurat vs singleR agreement in one figure ----
# built from the per-tissue CSVs written above; rows = tissues, columns = species
# panel order: rows of the combined figure
agree_all = map_dfr(c('colon', 'bone_marrow', 'fetal_liver', 'breast_cancer'), \(t)
  read_csv(file.path('figures', t, 'annotation', 'seurat_singleR_agreement.csv'), show_col_types = FALSE) |>
    mutate(tissue = t))
# breast cancer mouse is left blank (title only): its reference has no Cancer Epithelial type
p = purrr::pmap(distinct(agree_all, tissue, species), function(tissue, species) {
  title = sprintf('%s - %s', str_to_sentence(str_replace(tissue, '_', ' ')), species)
  if (tissue == 'breast_cancer' && species == 'mouse')
    return(ggplot() + labs(title = title) + theme_void(base_size = 9) + theme(plot.title = element_text(face = 'bold')))
  plot_annotation_agreement(filter(agree_all, tissue == !!tissue, species == !!species) |> select(-tissue),
                            title = title, base_size = 9, text_size = 2.3, axis_titles = c('SingleR (% of cells)', 'Seurat (% of cells)'))
})
p = patchwork::wrap_plots(p, ncol = 2, guides = 'collect') +
            patchwork::plot_annotation(
              tag_levels = 'A') &   # explanatory notes live in the thesis caption, not the figure
            theme(legend.position = 'bottom', plot.tag = element_text(face = 'bold', size = 13))
save_plot(p, 'figures/cross_tissue/seurat_singleR_agreement.png', width = 13, height = 18)
ggsave('figures/cross_tissue/seurat_singleR_agreement.pdf', p, width = 13, height = 18, device = cairo_pdf)   # vector, for the thesis

# heldout_auroc.png and heldout_auroc_celltype.png: see analysis/vote_vs_ebcc.R
