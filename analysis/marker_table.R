# Table of the marker genes for each cell type, tissue and species:  Rscript analysis/marker_table.R
# Marker lists are read from each tissue's annotation script exactly as in analysis/annotation_auroc_heldout.R (same files,
# variables and breast-cancer label maps), so these are the markers used for the marker AUROC (best annotator, EBCC
# checks). Panel membership comes from figures/cross_tissue/marker_coverage/<tissue>_<species>.csv
# (analysis/check_marker_coverage.R: markers specified but not in the dataset's gene panel); only markers on the panel can
# contribute to a score. Reads text and csv files only (no Seurat).
#   Rscript --vanilla analysis/marker_table.R
# Output: figures/cross_tissue/marker_table.csv (one row per tissue, cell type, species: on panel, not on panel),
#         figures/cross_tissue/marker_table.tex (longtable: tissue, cell type, human markers, mouse markers; on panel)
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(stringr)})

# ---------------- config ----------------
cfg = list(   # as analysis/annotation_auroc_heldout.R
  bone_marrow = list(human = list(markers = 'bone_marrow/02_human_annotate.R'), mouse = list(markers = 'bone_marrow/02_mouse_annotate.R')),
  colon = list(human = list(markers = 'colon/02_human_annotated.R'), mouse = list(markers = 'colon/02_mouse_annotate.qmd')),
  breast_cancer = list(
    human = list(markers = 'breast_cancer/02_human_annotate.qmd',
                 label_map = c('B-cells' = 'B', 'CAFs' = 'CAF', 'Cancer Epithelial' = 'Tumour', 'Endothelial' = 'Endothelial',
                               'Myeloid' = 'Myeloid', 'Normal Epithelial' = 'Basal', 'Plasmablasts' = 'Plasmablast',
                               'PVL' = 'PVL', 'T-cells' = 'T_NK')),
    mouse = list(markers = 'breast_cancer/02_mouse_annotate.qmd', markers_var = 'markers_pymt',
                 label_map = c('B-cells' = 'B_cells', 'CAFs' = 'CAF', 'Cancer Epithelial' = 'Tumour', 'Endothelial' = 'Endothelial',
                               'Myeloid' = 'Myeloid', 'Normal Epithelial' = 'Basal', 'Plasmablasts' = 'Plasma',
                               'PVL' = 'PVL', 'T-cells' = 'T_cells'))),
  fetal_liver = list(human = list(markers = 'fetal_liver/02_human_annotate.R'), mouse = list(markers = 'fetal_liver/02_mouse_annotate.R')))
tissue_labels = c(colon = 'Colon', bone_marrow = 'Bone marrow', fetal_liver = 'Fetal liver', breast_cancer = 'Breast cancer')
coverage_dir = 'figures/cross_tissue/marker_coverage'
out = 'figures/cross_tissue'
# -----------------------------------------

read_markers = function(script, var = 'markers') {   # as analysis/annotation_auroc_heldout.R / analysis/visualisations.R
  l = readLines(script, warn = FALSE)
  s = grep(paste0('^', var, ' = list'), l)[1]
  e = s - 1 + which(l[s:length(l)] == ')')[1]
  env = new.env()
  eval(parse(text = l[s:e]), envir = env)
  env[[var]]
}
annotation_markers = function(a) {
  m = read_markers(a$markers, if (is.null(a$markers_var)) 'markers' else a$markers_var)
  if (!is.null(a$label_map)) m = setNames(m[a$label_map], names(a$label_map))
  m[!vapply(m, is.null, logical(1))]
}

tab = bind_rows(lapply(names(cfg), \(t) bind_rows(lapply(names(cfg[[t]]), function(s) {
  m = annotation_markers(cfg[[t]][[s]])
  cov = read_csv(file.path(coverage_dir, sprintf('%s_%s.csv', t, s)), show_col_types = FALSE)
  stopifnot(setequal(names(m), cov$type))   # same cell types as the coverage check
  bind_rows(lapply(names(m), function(ty) {
    missing = str_split(cov$missing[cov$type == ty], ';')[[1]]
    missing = missing[!is.na(missing) & missing != '']
    stopifnot(all(missing %in% m[[ty]]), length(m[[ty]]) == cov$n_specified[cov$type == ty])
    tibble(tissue = t, species = s, type = ty, n_specified = length(m[[ty]]),
           on_panel = paste(setdiff(m[[ty]], missing), collapse = ', '), n_on_panel = length(setdiff(m[[ty]], missing)),
           not_on_panel = paste(missing, collapse = ', '))
  }))
}))))
# cell types are matched between species by a normalised name, as in analysis/annotation_auroc_heldout.R (key()); shown with
# the human (EBCC) name, or the mouse name with '_' as spaces for mouse-only types
key = function(x) tolower(gsub('[^A-Za-z0-9]+', '_', x))
tab = tab |> mutate(k = key(type)) |> group_by(tissue, k) |>
  mutate(label = if (any(species == 'human')) type[species == 'human'][1] else gsub('_', ' ', type[1])) |> ungroup() |>
  mutate(tissue = factor(tissue, levels = names(tissue_labels))) |> arrange(tissue, label, species) |> select(-k) |> relocate(label, .after = tissue)
write_csv(tab, file.path(out, 'marker_table.csv'))

# LaTeX longtable: markers on the panel, human and mouse side by side
esc = \(x) str_replace_all(x, c('_' = '\\\\_', '&' = '\\\\&', '%' = '\\\\%', '#' = '\\\\#'))
w = tab |> select(tissue, label, species, on_panel, n_on_panel, n_specified) |>
  mutate(cell = ifelse(n_on_panel == 0, '--', esc(on_panel))) |> select(tissue, label, species, cell) |>
  pivot_wider(names_from = species, values_from = cell, values_fill = '--')
lines = c('{\\footnotesize',
          '\\begin{longtable}{@{} p{0.12\\textwidth} p{0.16\\textwidth} p{0.33\\textwidth} p{0.33\\textwidth} @{}}',
          '\\caption{Marker genes used for each cell type, by tissue and species (on the gene panel of each dataset; markers specified in the annotation but absent from the panel are omitted, and -- marks a cell type with no measured marker or no marker list in that species).}',
          '\\label{tab:markers} \\\\',
          '\\toprule', 'Tissue & Cell type & Human markers & Mouse markers \\\\', '\\midrule', '\\endfirsthead',
          '\\toprule', 'Tissue & Cell type & Human markers & Mouse markers \\\\', '\\midrule', '\\endhead',
          '\\midrule \\multicolumn{4}{r}{\\textit{continued on next page}} \\\\', '\\endfoot',
          '\\bottomrule', '\\endlastfoot')
for (t in levels(w$tissue)) {
  x = filter(w, tissue == t)
  for (i in seq_len(nrow(x)))
    lines = c(lines, sprintf('%s & %s & %s & %s \\\\', if (i == 1) tissue_labels[[t]] else '', esc(x$label[i]), x$human[i], x$mouse[i]))
  if (t != tail(levels(w$tissue), 1)) lines = c(lines, '\\midrule')
}
lines = c(lines, '\\end{longtable}', '}')
writeLines(lines, file.path(out, 'marker_table.tex'))

s = tab |> group_by(tissue, species) |> summarise(types = n(), specified = sum(n_specified), on_panel = sum(n_on_panel), .groups = 'drop')
cat('Marker genes per tissue and species (specified in the annotation script; on the gene panel):\n')
print(as.data.frame(s), row.names = FALSE)
cat('wrote', file.path(out, 'marker_table.{csv,tex}'), '\n')
