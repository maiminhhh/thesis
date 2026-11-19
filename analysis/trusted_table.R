# LaTeX table of the trusted cell types of each annotator, per tissue and species, with the trust weight frac used:
#   Rscript analysis/trusted_table.R
# The `trusted` lists are read from each species' annotation script (the list passed to the EBCC prior: diagonal
# 4 + frac * n_k / M for each trusted (annotator, class); EBCC_PIPELINE.md section 2); frac is the value finally used
# (data/02_annotated_data/<species>_ebcc_frac_auroc_chosen.csv). With frac = 0 the lists have no effect.
# Rows: every EBCC label of the tissue (both species) plus any trusted class that is not an EBCC label.
#   Rscript --vanilla analysis/trusted_table.R
# Output: figures/cross_tissue/trusted_table.{tex,csv}
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(stringr)})

# ---------------- config ----------------
scripts = list(
  colon         = c(human = 'colon/02_human_annotated.R',         mouse = 'colon/02_mouse_annotate.qmd'),
  bone_marrow   = c(human = 'bone_marrow/02_human_annotate.R',    mouse = 'bone_marrow/02_mouse_annotate.R'),
  fetal_liver   = c(human = 'fetal_liver/02_human_annotate.R',    mouse = 'fetal_liver/02_mouse_annotate.R'),
  breast_cancer = c(human = 'breast_cancer/02_human_annotate.qmd', mouse = 'breast_cancer/02_mouse_annotate.qmd'))
tissue_labels = c(colon = 'Colon', bone_marrow = 'Bone marrow', fetal_liver = 'Fetal liver', breast_cancer = 'Breast cancer')
annotators = c(celltype_seurat = 'S', celltype_singleR = 'R', celltype_seurat_cross = 'S$_\\times$', celltype_singleR_cross = 'R$_\\times$')
labels_file = 'figures/dissection/dissection_cells.csv.gz'   # EBCC labels of every cell (analysis/dissection_figure_compute.py)
out = 'figures/cross_tissue'
# -----------------------------------------

read_list = function(script, var = 'trusted') {   # as read_markers in analysis/annotation_auroc_heldout.R
  l = readLines(script, warn = FALSE)
  s = grep(paste0('^', var, ' *= *list'), l)[1]
  e = s - 1 + which(l[s:length(l)] == ')')[1]
  env = new.env(); eval(parse(text = l[s:e]), envir = env); env[[var]]
}
tr = bind_rows(lapply(names(scripts), \(t) bind_rows(lapply(names(scripts[[t]]), function(s) {
  x = read_list(scripts[[t]][[s]])
  stopifnot(setequal(names(x), names(annotators)))
  bind_rows(lapply(names(x), \(a) tibble(tissue = t, species = s, annotator = a, label = x[[a]])))
}))))
frac = bind_rows(lapply(names(scripts), \(t) bind_rows(lapply(c('human', 'mouse'), \(s)
  read_csv(sprintf('%s/data/02_annotated_data/%s_ebcc_frac_auroc_chosen.csv', t, s), show_col_types = FALSE) |>
    transmute(tissue = t, species = s, M, frac)))))
ebcc = read_csv(labels_file, show_col_types = FALSE, col_select = c(tissue, label)) |> distinct()

rows = bind_rows(ebcc, distinct(tr, tissue, label)) |> distinct() |> mutate(ebcc = paste(tissue, label) %in% paste(ebcc$tissue, ebcc$label))
grid = rows |> crossing(species = c('human', 'mouse'), annotator = names(annotators)) |>
  left_join(tr |> mutate(trusted = TRUE), by = c('tissue', 'label', 'species', 'annotator')) |>
  mutate(trusted = !is.na(trusted), tissue = factor(tissue, levels = names(tissue_labels))) |>
  arrange(tissue, label, species, match(annotator, names(annotators)))
write_csv(grid |> left_join(frac, by = c('tissue', 'species')), file.path(out, 'trusted_table.csv'))

esc = \(x) str_replace_all(x, c('_' = '\\\\_', '&' = '\\\\&', '%' = '\\\\%'))
fmt_frac = \(v) ifelse(v == round(v), sprintf('%d', as.integer(v)), sprintf('%g', v))
hdr = c(sprintf('\\multicolumn{2}{c}{} & \\multicolumn{4}{c}{Human} & \\multicolumn{4}{c}{Mouse} \\\\'),
        '\\cmidrule(lr){3-6} \\cmidrule(l){7-10}',
        sprintf('Tissue & Cell type & %s & %s \\\\', paste(annotators, collapse = ' & '), paste(annotators, collapse = ' & ')))
lines = c('{\\footnotesize',
          '\\begin{longtable}{@{} l l cccc cccc @{}}',
          paste0('\\caption{Cell types trusted for each annotator in the EBCC prior, per tissue and species. ',
                 'S: Seurat and R: SingleR with a same-species reference; S$_\\times$ and R$_\\times$: with the other species\' reference. ',
                 'The prior adds $\\mathit{frac}\\cdot n_k/M$ to the diagonal for each trusted class $k$ of an annotator ($n_k$: cells it gave that label; $M$: subtypes per class); ',
                 'with $\\mathit{frac}=0$ the lists have no effect.',
                 if (any(grid$trusted & !grid$ebcc)) ' Cell types in italics are not EBCC labels.' else '', '}'),
          '\\label{tab:trusted} \\\\', '\\toprule', hdr, '\\midrule', '\\endfirsthead',
          '\\toprule', hdr, '\\midrule', '\\endhead',
          '\\midrule \\multicolumn{10}{r}{\\textit{continued on next page}} \\\\', '\\endfoot',
          '\\bottomrule', '\\endlastfoot')
for (t in levels(grid$tissue)) {
  g = filter(grid, tissue == t)
  labs = g |> distinct(label, ebcc) |> arrange(!ebcc, label)
  for (i in seq_len(nrow(labs))) {
    cells = sapply(c('human', 'mouse'), \(s) sapply(names(annotators), \(a)
      if (filter(g, label == labs$label[i], species == s, annotator == a)$trusted) '\\checkmark' else ''))
    nm = if (labs$ebcc[i]) esc(labs$label[i]) else sprintf('\\textit{%s}', esc(labs$label[i]))
    lines = c(lines, sprintf('%s & %s & %s \\\\', if (i == 1) tissue_labels[[t]] else '', nm, paste(cells, collapse = ' & ')))
  }
  f = filter(frac, tissue == t)
  lines = c(lines, sprintf(' & $\\mathit{frac}$ & \\multicolumn{4}{c}{%s} & \\multicolumn{4}{c}{%s} \\\\',
                           fmt_frac(f$frac[f$species == 'human']), fmt_frac(f$frac[f$species == 'mouse'])))
  if (t != tail(levels(grid$tissue), 1)) lines = c(lines, '\\midrule')
}
lines = c(lines, '\\end{longtable}', '}')
writeLines(lines, file.path(out, 'trusted_table.tex'))
print(as.data.frame(frac), row.names = FALSE)
cat('trusted (annotator, class) pairs per tissue and species:\n'); print(as.data.frame(count(tr, tissue, species)), row.names = FALSE)
cat('trusted classes that are not EBCC labels:\n'); print(as.data.frame(distinct(filter(grid, trusted, !ebcc), tissue, label)), row.names = FALSE)
cat('wrote', file.path(out, 'trusted_table.{tex,csv}'), '\n')
