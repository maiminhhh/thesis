# Marker AUROC of every labelling on the EBCC evaluation half:  Rscript analysis/annotation_auroc_heldout.R
# Labellings: the four annotators, majority vote, EBCC with frac = 0 and EBCC with the chosen frac.
# Evaluation half: each tissue's 02 frac selection split the cells of the 02 annotated object 50/50 with
# set.seed(1); sample(n, floor(n / 2)) ('selection'); the other half is reported. It is redrawn here the
# same way, and the final EBCC labels must reproduce the stored chosen-frac AUROC on it (checked).
# Marker score and AUROC as in analysis/visualisations.R: mean log-normalised expression of a type's panel
# markers; AUROC per type over the cells an annotator labelled (SingleR leaves some unlabelled);
# summary = mean over the types the labelling calls. EBCC frac = 0 is read from the stored selection
# tables (its fit is not kept for every tissue); for chosen frac the stored value is shown, after the check.
# Majority vote: most frequent of the four labels a cell has; ties broken at random (seed 1).
# Output: figures/cross_tissue/annotation_auroc_heldout.{csv,png}
source('library.R')

annotation = list(
  bone_marrow = list(
    human = list(obj = 'bone_marrow/data/02_annotated_data/02_human_annotated.rds', markers = 'bone_marrow/02_human_annotate.R'),
    mouse = list(obj = 'bone_marrow/data/02_annotated_data/02_mouse_annotated.rds', markers = 'bone_marrow/02_mouse_annotate.R')),
  colon = list(
    human = list(obj = 'colon/data/02_annotated_data/02_human_annotated.rds', markers = 'colon/02_human_annotated.R'),
    mouse = list(obj = 'colon/data/02_annotated_data/02_mouse_annotated.rds', markers = 'colon/02_mouse_annotate.qmd')),
  breast_cancer = list(
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
                               'PVL' = 'PVL', 'T-cells' = 'T_cells'))),
  fetal_liver = list(
    human = list(obj = 'fetal_liver/data/02_annotated_data/02_human_annotated.rds', markers = 'fetal_liver/02_human_annotate.R'),
    mouse = list(obj = 'fetal_liver/data/02_annotated_data/02_mouse_annotated.rds', markers = 'fetal_liver/02_mouse_annotate.R')))
annotators = c(Seurat = 'celltype_seurat', SingleR = 'celltype_singleR',
               'Seurat (cross)' = 'celltype_seurat_cross', 'SingleR (cross)' = 'celltype_singleR_cross')

# same helpers as analysis/visualisations.R
read_markers = function(script, var = 'markers') {
  l = readLines(script)
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
auroc = function(score, pos) {            # Mann-Whitney form, ties averaged
  r = rank(score); n1 = as.numeric(sum(pos)); n0 = as.numeric(sum(!pos))
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
key = function(x) tolower(gsub('[^A-Za-z0-9]+', '_', x))

# mean AUROC over the types a labelling calls, on the cells it labelled
mean_auroc = function(lab, scores) {
  lab = key(as.character(lab)); ok = !is.na(lab)
  a = vapply(colnames(scores), function(k) {
    pos = lab[ok] == key(k)
    if (!any(pos) || all(pos)) NA_real_ else auroc(scores[ok, k], pos)
  }, numeric(1))
  c(mean = mean(a, na.rm = TRUE), n_types = sum(!is.na(a)), n_cells = sum(ok))
}

majority_vote = function(meta) {
  L = sapply(annotators, function(col) as.character(meta[[col]]))
  classes = sort(unique(na.omit(as.vector(L))))
  counts = sapply(classes, function(k) rowSums(L == k, na.rm = TRUE))
  set.seed(1)
  out = classes[max.col(counts, ties.method = 'random')]
  out[rowSums(counts) == 0] = NA
  n_max = apply(counts, 1, max)
  attr(out, 'tied') = mean(rowSums(counts == n_max) > 1 & n_max > 0)
  out
}

stored = function(tissue, species) {
  d = read.csv(sprintf('%s/data/02_annotated_data/%s_ebcc_frac_auroc.csv', tissue, species))
  held = if ('mean_auroc_report' %in% names(d)) d$mean_auroc_report else d$mean_auroc_heldout
  list(f0 = held[d$frac == 0], chosen = held[d$chosen %in% c(TRUE, 'TRUE')],
       frac = d$frac[d$chosen %in% c(TRUE, 'TRUE')])
}

res = list(); checks = list()
for (tissue in names(annotation)) for (species in names(annotation[[tissue]])) {
  a = annotation[[tissue]][[species]]
  message(sprintf('%s %s: loading', tissue, species))
  obj = readRDS(a$obj)
  assay = DefaultAssay(obj)
  rownames(obj) = obj[[assay]]@meta.data$gene_name
  X = LayerData(obj, assay = assay, layer = 'data')
  meta = obj@meta.data
  rm(obj); gc()
  markers = lapply(annotation_markers(a), intersect, rownames(X))
  markers = markers[lengths(markers) > 0]
  scores = sapply(markers, function(g) Matrix::colMeans(X[g, , drop = FALSE]))
  rm(X); gc()

  set.seed(1)                                     # the 02 frac-selection split
  sel = seq_len(nrow(meta)) %in% sample(nrow(meta), floor(nrow(meta) / 2))
  ev = !sel

  mv = majority_vote(meta)
  labellings = c(lapply(annotators, function(col) meta[[col]]),
                 list('Majority vote' = mv, 'EBCC (check)' = meta$ebcc_label))
  r = imap_dfr(labellings, function(lab, name) {
    m = mean_auroc(lab[ev], scores[ev, , drop = FALSE])
    tibble(labelling = name, mean_auroc = m[['mean']], n_types = m[['n_types']], n_cells = m[['n_cells']])
  })
  st = stored(tissue, species)
  check = r$mean_auroc[r$labelling == 'EBCC (check)']
  checks[[length(checks) + 1]] = tibble(tissue, species, recomputed = check, stored = st$chosen,
                                        diff = check - st$chosen, mv_tied = attr(mv, 'tied'))
  r = bind_rows(filter(r, labelling != 'EBCC (check)'),
                tibble(labelling = c('EBCC (f = 0)', 'EBCC (chosen f)'), mean_auroc = c(st$f0, st$chosen),
                       n_types = NA, n_cells = sum(ev)))
  res[[length(res) + 1]] = mutate(r, tissue = tissue, species = species, chosen_f = st$frac, .before = 1)
  message(sprintf('  EBCC check: recomputed %.4f vs stored %.4f; majority-vote ties %.1f%%',
                  check, st$chosen, 100 * attr(mv, 'tied')))
}
res = bind_rows(res); checks = bind_rows(checks)
print(checks)
stopifnot(all(abs(checks$diff) < 1e-6))            # the evaluation half is the stored one

dir.create('figures/cross_tissue', showWarnings = FALSE, recursive = TRUE)
write_csv(res, 'figures/cross_tissue/annotation_auroc_heldout.csv')
write_csv(checks, 'figures/cross_tissue/annotation_auroc_heldout_checks.csv')

tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
f_lab = res |> distinct(tissue, species, chosen_f) |>
  pivot_wider(names_from = species, values_from = chosen_f) |>
  mutate(strip = sprintf('%s\n(chosen f: human %s, mouse %s)', tissue_labels[tissue], human, mouse))
d = res |> left_join(select(f_lab, tissue, strip), by = 'tissue') |>
  mutate(strip = factor(strip, levels = f_lab$strip[match(names(tissue_labels), f_lab$tissue)]),
         labelling = factor(labelling, levels = c(names(annotators), 'Majority vote', 'EBCC (f = 0)', 'EBCC (chosen f)')),
         kind = if_else(str_starts(as.character(labelling), 'EBCC'), 'EBCC', 'other'),
         species = factor(species, levels = c('human', 'mouse')))
p = ggplot(d, aes(labelling, mean_auroc, colour = species, shape = species)) +
  geom_rect(data = distinct(d, strip), aes(xmin = 5.5, xmax = 7.5, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = 'grey93') +
  geom_point(size = 3, position = position_dodge(width = 0.35)) +
  facet_wrap(~ strip, nrow = 1) +
  scale_colour_manual(values = c(human = '#2a78d6', mouse = '#eb6834'), name = NULL) +
  scale_shape_manual(values = c(human = 16, mouse = 17), name = NULL) +
  labs(x = NULL, y = 'Mean marker AUROC (evaluation half)',
       title = 'Agreement of each labelling with panel markers',
       subtitle = 'Evaluation half of the cells (not used to choose f). Shaded: EBCC.',
       caption = paste('Mean over the cell types each labelling calls; AUROC of the mean log-normalised expression',
                       'of the type\'s panel markers. Majority vote: ties broken at random.')) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), legend.position = 'top', plot.title.position = 'plot',
        strip.background = element_rect(fill = 'grey95'), plot.caption = element_text(hjust = 0, colour = 'grey35'))
ggsave('figures/cross_tissue/annotation_auroc_heldout.png', p, width = 14, height = 5.8, dpi = 300, bg = 'white')
print(res |> select(tissue, species, labelling, mean_auroc) |> mutate(mean_auroc = round(mean_auroc, 3)) |>
        pivot_wider(names_from = labelling, values_from = mean_auroc), width = Inf)
