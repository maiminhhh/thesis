# Majority vote vs EBCC (and the single annotators) scored fairly. Run from the repository root.
# The held-out AUROC figure (analysis/visualisations.R, heldout_auroc.png) scores each labelling
# on the cells it labels, so the vote skips its ties (the hardest cells) while EBCC
# labels every cell. Two comparisons, both on the evaluation half (EBCC section 3b
# split, set.seed(1)):
#   common: cells every one of the 7 labellings labels (vote not tied, SingleR not
#           pruned); fair, but on the easier cells
#   all:    every evaluation cell, for the labellings that label all of them (Seurat,
#           Seurat cross, EBCC f = 0, EBCC chosen f); includes the hard cells
# Per labelling: mean marker AUROC over the types it calls in that cell set (as in
# heldout_auroc.png). 95% intervals from B bootstrap resamples of the cells; the
# difference to EBCC (chosen f) is paired (same resample).
#
# Usage:  Rscript analysis/vote_vs_ebcc.R <tissue> <species>   one dataset -> figures/cross_tissue/vote_vs_ebcc/<tissue>_<species>.csv
#         Rscript analysis/vote_vs_ebcc.R                      combine the CSVs -> figures/cross_tissue/vote_vs_ebcc*.png/.pdf/.csv,
#                                                     heldout_auroc.png, heldout_auroc_celltype.png (+ .pdf)
#         VOTE_B=0 Rscript analysis/vote_vs_ebcc.R <tissue> <species>   per-cell-type CSV only (skips the bootstrap)
# bone_marrow's f = 0 fit was not cached by section 3b; refit it first with bone_marrow/ebcc_refit_f0.R (same call, seed 1)
# into bone_marrow/data/02_annotated_data/ebcc_frac_fits/<species>_frac0.rds.

# tissue config and helpers (tissues, annotators, annotation_markers, majority_vote, label_sets,
# save_plot, species_cols): every top-level expression of analysis/visualisations.R before its run loop
ex = parse('analysis/visualisations.R')
run_at = which(startsWith(sapply(ex, \(e) deparse(e)[1]), 'for (tissue in names(tissues))'))
# (skipping its run_only filter, which reads this script's <tissue> <species> arguments as tissue names)
for (e in ex[seq_len(run_at - 1)]) {
  if (grepl('run_only', paste(deparse(e), collapse = ' '))) next
  eval(e, globalenv())
}

B = as.integer(Sys.getenv('VOTE_B', '200'))   # VOTE_B=0: per-cell-type AUROC only (no bootstrap)
tissue_order = c('colon', 'bone_marrow', 'fetal_liver', 'breast_cancer')
out_dir = 'figures/cross_tissue/vote_vs_ebcc'
key = function(x) tolower(gsub('[^A-Za-z0-9]+', '_', x))
all_cell_sets = c('Seurat', 'Seurat (cross)', 'EBCC (f = 0)', 'EBCC (chosen f)')

# cached EBCC fit for (tissue, species, frac)
fit_file = function(tissue, species, f, M) {
  d = file.path(tissue, 'data/02_annotated_data')
  switch(tissue,
         bone_marrow   = file.path(d, 'ebcc_frac_fits', sprintf('%s_frac%s.rds', species, format(f))),
         colon         = file.path(d, 'ebcc_frac_fits', sprintf('%s_frac%s.rds', species, format(f))),
         breast_cancer = file.path(d, sprintf('%s_ebcc_frac_fits', species), sprintf('M%s_frac%s.rds', M, format(f))),
         fetal_liver   = file.path(d, 'ebcc_frac', species, sprintf('frac%s.rds', format(f))))
}
fit_labels = function(file, cells) {
  x = readRDS(file)
  lab = as.character(x$label)
  names(lab) = if (!is.null(names(x$label))) names(x$label) else rownames(x$posterior)
  stopifnot(all(cells %in% names(lab)))
  unname(lab[cells])
}

# AUROC from precomputed ranks (Mann-Whitney, ties averaged by rank())
auc_rank = function(r, pos) {
  n1 = as.numeric(sum(pos)); n0 = as.numeric(length(pos) - n1)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# per labelling: mean AUROC over its called types (fixed on the full cell set), point
# estimate + bootstrap; idx = cell rows of `score`, labs = keyed labels on those rows
score_set = function(score, labs, idx, set_name) {
  types = lapply(labs, \(l) colnames(score)[vapply(colnames(score), \(k)
    { p = l[idx] == key(k); any(p) && !all(p) }, logical(1))])
  means = function(i) {
    r = lapply(setNames(nm = unique(unlist(types))), \(k) rank(score[i, k]))   # one rank per type
    vapply(names(labs), \(n) mean(vapply(types[[n]], \(k) auc_rank(r[[k]], labs[[n]][i] == key(k)),
                                         numeric(1)), na.rm = TRUE), numeric(1))
  }
  est = means(idx)
  boot = t(replicate(B, means(sample(idx, replace = TRUE))))
  ref = boot[, 'EBCC (chosen f)']
  ci = function(x) quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  tibble(set = set_name, labelling = names(labs), n_cells = length(idx), n_types = lengths(types),
         mean_auroc = est, lo = apply(boot, 2, \(x) ci(x)[1]), hi = apply(boot, 2, \(x) ci(x)[2]),
         diff_vs_ebcc = est - est[['EBCC (chosen f)']],
         diff_lo = apply(boot - ref, 2, \(x) ci(x)[1]), diff_hi = apply(boot - ref, 2, \(x) ci(x)[2]))
}

compare_one = function(tissue, species) {
  t0 = Sys.time()
  tick = function(msg) message(sprintf('[%s %s %5.1f min] %s', tissue, species,
                                       as.numeric(Sys.time() - t0, units = 'mins'), msg))
  a = tissues[[tissue]]$annotation[[species]]
  obj = readRDS(a$obj)
  assay = DefaultAssay(obj)
  rownames(obj) = obj[[assay]]@meta.data$gene_name
  X = LayerData(obj, assay = assay, layer = 'data')
  meta = obj@meta.data
  cells = colnames(obj)
  rm(obj); gc(); tick('object loaded')
  set.seed(1)
  held = !seq_len(ncol(X)) %in% sample(ncol(X), floor(ncol(X) / 2))
  all_markers = lapply(annotation_markers(a), intersect, rownames(X))
  markers = all_markers[lengths(all_markers) > 0]
  score = sapply(markers, \(g) Matrix::colMeans(X[g, held, drop = FALSE]))   # held cells x types
  rm(X); gc()

  ch = read_csv(file.path(dirname(a$obj), str_glue('{species}_ebcc_frac_auroc_chosen.csv')), show_col_types = FALSE)
  # bone_marrow: the chosen fit (f = 5) is the object's ebcc_label; others from the cache
  chosen = if (tissue == 'bone_marrow') as.character(meta$ebcc_label) else
    fit_labels(fit_file(tissue, species, ch$frac, ch$M), cells)
  chosen_match = mean(chosen == as.character(meta$ebcc_label))   # cached chosen fit vs final labels (expect 1)
  raw = list('Seurat'          = meta[[annotators[['seurat']]]],
              'SingleR'         = meta[[annotators[['singleR']]]],
              'Seurat (cross)'  = meta[[annotators[['seurat_cross']]]],
              'SingleR (cross)' = meta[[annotators[['singleR_cross']]]],
              'Majority vote'   = majority_vote(meta),
              'EBCC (f = 0)'    = fit_labels(fit_file(tissue, species, 0, ch$M), cells),
              'EBCC (chosen f)' = chosen) |>
    lapply(as.character)
  labs = lapply(raw, \(l) key(l[held]))
  tick('labels ready')

  labelled = sapply(labs, \(l) !is.na(l))                       # held cells x labellings
  common = which(rowSums(!labelled) == 0)
  all_idx = which(rowSums(!labelled[, all_cell_sets]) == 0)     # = all held cells unless Seurat has NAs
  info = \(d) mutate(d, tissue = tissue, species = species, chosen_frac = ch$frac, n_held = sum(held),
                     .before = 1)

  # per cell type on the common cells, for the heatmap: every marker group (n_markers = on
  # the panel; 0 -> cannot be scored) plus labels without one
  extra = unique(na.omit(unlist(raw)))
  extra = extra[!key(extra) %in% key(names(all_markers))]
  types = tibble(celltype = c(names(all_markers), extra),
                 n_markers = c(lengths(all_markers), rep(0L, length(extra))))
  per_class = imap_dfr(labs, function(l, set) {
    map_dfr(seq_len(nrow(types)), function(i) {
      k = types$celltype[i]; pos = l[common] == key(k)
      tibble(label_set = set, celltype = k, n_markers = types$n_markers[i], n_cells = sum(pos),
             auroc = if (types$n_markers[i] == 0 || !any(pos) || all(pos)) NA_real_ else
               auroc(score[common, k], pos))
    })
  }) |> info() |> mutate(n_common = length(common))

  res = NULL
  if (B > 0) {
    set.seed(2)
    res = bind_rows(score_set(score, labs, common, 'common'),
                    score_set(score, labs[all_cell_sets], all_idx, 'all')) |>
      info() |> mutate(chosen_match = chosen_match)
    tick('bootstrap done')
  }
  list(summary = res, per_class = per_class)
}

# mean AUROC per labelling with 95% bootstrap intervals, one facet per tissue; the labelling used
# (EBCC, chosen f) is shaded, bold and labelled with its values, its dashed line a reference
plot_auroc = function(d, levels, title = NULL, xlab = 'Mean marker AUROC, same cells per panel') {
  chosen = 'EBCC (chosen f)'
  d = d |> mutate(tissue = factor(tissue, levels = tissue_order), labelling = factor(labelling, levels = rev(levels)))
  ref = filter(d, labelling == chosen)
  strip = ref |> distinct(tissue, species, chosen_frac, n_cells, n_held) |> group_by(tissue) |>
    summarise(lab = sprintf('%s\n%s', str_to_sentence(str_replace(tissue[1], '_', ' ')),
                            paste(sprintf('%s: %s cells (%.0f%%), f = %g', species, scales::comma(n_cells),
                                          100 * n_cells / n_held, chosen_frac), collapse = '\n')))
  row = which(rev(levels) == chosen)
  dodge = position_dodge(width = 0.6)
  ggplot(d, aes(mean_auroc, labelling, colour = species)) +
    geom_rect(data = data.frame(ymin = row - 0.5, ymax = row + 0.5), inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = 'EBCC (chosen f): labels used in the benchmark')) +
    geom_vline(data = ref, aes(xintercept = mean_auroc, colour = species, linetype = 'EBCC (chosen f), per species'),
               linewidth = 0.35, alpha = 0.7, show.legend = c(colour = FALSE, linetype = TRUE)) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = dodge, size = 0.2, linewidth = 0.7) +
    geom_text(data = ref, aes(x = hi, label = sprintf('%.3f', mean_auroc)),
              position = dodge, hjust = -0.25, size = 2.7, fontface = 'bold', show.legend = FALSE) +
    facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(setNames(strip$lab, strip$tissue))) +
    scale_colour_manual(values = c(human = species_cols[['hsapiens']], mouse = species_cols[['mmusculus']]),
                        name = sprintf('Mean, 95%% bootstrap CI\n(%d resamples of cells)', B)) +
    scale_linetype_manual(values = 'dashed', name = 'Reference line') +
    scale_fill_manual(values = '#FCEFC7', name = 'Shaded row') +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2, override.aes = list(colour = 'grey30', alpha = 1)),
           fill = guide_legend(order = 3)) +
    scale_x_continuous(breaks = seq(0, 1, 0.1), expand = expansion(mult = c(0.05, 0.18))) +
    labs(title = title, x = xlab, y = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'top', legend.title = element_text(size = 10), legend.key.width = unit(1.2, 'lines'),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          axis.text.y = element_text(face = ifelse(rev(levels) == chosen, 'bold', 'plain')),
          strip.background = element_blank(), strip.text = element_text(hjust = 0, size = 9),
          plot.title = element_text(face = 'bold'), plot.caption = element_text(hjust = 0, size = 9))
}

args = commandArgs(TRUE)
if (length(args) == 2) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  r = compare_one(args[1], args[2])
  if (!is.null(r$summary)) write_csv(r$summary, file.path(out_dir, sprintf('%s_%s.csv', args[1], args[2])))
  write_csv(r$per_class, file.path(out_dir, sprintf('%s_%s_per_class.csv', args[1], args[2])))
  message('done ', args[1], ' ', args[2])
  quit(save = 'no')
}

# ---- combine: one CSV, one figure per comparison ----
res = map_dfr(tissue_order, \(t) map_dfr(c('human', 'mouse'), \(sp)
  read_csv(file.path(out_dir, sprintf('%s_%s.csv', t, sp)), show_col_types = FALSE)))
write_csv(res, 'figures/cross_tissue/vote_vs_ebcc.csv')

cap = 'Evaluation half of cells (not used to choose f); strips: cells used (% of the half) and chosen f. Mean AUROC over the types each labelling calls.'
p1 = plot_auroc(filter(res, set == 'common'), label_sets,
                'Common cells: labelled by all seven labellings (vote not tied, SingleR not pruned)') + labs(caption = cap)
p2 = plot_auroc(filter(res, set == 'all'), all_cell_sets, 'All cells: labellings that label every cell') + labs(caption = cap)
save_plot(p1, 'figures/cross_tissue/vote_vs_ebcc_common.png', width = 14, height = 5.5)
ggsave('figures/cross_tissue/vote_vs_ebcc_common.pdf', p1, width = 14, height = 5.5)
save_plot(p2, 'figures/cross_tissue/vote_vs_ebcc_all.png', width = 14, height = 4)
ggsave('figures/cross_tissue/vote_vs_ebcc_all.pdf', p2, width = 14, height = 4)

# ---- heldout_auroc.png: mean AUROC per labelling on the common cells ----
p = plot_auroc(filter(res, set == 'common'), label_sets, xlab = 'Mean marker AUROC (common cells of the evaluation half)') +
  labs(caption = paste('Same cells for every labelling: evaluation half (not used to choose f), restricted to cells all seven label (vote not tied, SingleR not pruned).',
                       'Strips: cells used (% of the half) and chosen f. Mean AUROC over the types each labelling calls.', sep = '\n'))
save_plot(p, 'figures/cross_tissue/heldout_auroc.png', width = 13, height = 5.5)
ggsave('figures/cross_tissue/heldout_auroc.pdf', p, width = 13, height = 5.5)

# ---- heldout_auroc_celltype.png: AUROC per cell type on the common cells ----
pc = map_dfr(tissue_order, \(t) map_dfr(c('human', 'mouse'), \(sp)
  read_csv(file.path(out_dir, sprintf('%s_%s_per_class.csv', t, sp)), show_col_types = FALSE)))
write_csv(pc, 'figures/cross_tissue/heldout_auroc_celltype.csv')
p = purrr::map(tissue_order, function(t) {
  r = filter(pc, tissue == t)
  cov = r |> distinct(species, n_common, n_held)
  plot_heldout_celltype(r, title = sprintf('%s  (%s)', str_to_sentence(str_replace(t, '_', ' ')),
                                           paste(sprintf('%s: %.0f%% of cells', cov$species,
                                                         100 * cov$n_common / cov$n_held), collapse = ', ')))
})
heights = purrr::map_int(tissue_order, \(t) n_distinct(key(pc$celltype[pc$tissue == t & pc$n_cells > 0])) + 1L)
p = patchwork::wrap_plots(p, ncol = 1, heights = heights, guides = 'collect') +
  patchwork::plot_annotation(
    tag_levels = 'A',
    caption = paste('Marker AUROC per cell type on the common cells of the evaluation half (labelled by all seven labellings; % per species in each title);',
                    'bottom row: mean over the types each labelling calls. Dark grey: no markers for the type on that species\' panel.',
                    '"–": labelling never calls the type on these cells. Blank: not a cell type in that species. EBCC f chosen on the other half.',
                    sep = '\n'),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8))) &
  theme(legend.position = 'bottom', plot.tag = element_text(face = 'bold', size = 13))
save_plot(p, 'figures/cross_tissue/heldout_auroc_celltype.png', width = 13, height = 22)
ggsave('figures/cross_tissue/heldout_auroc_celltype.pdf', p, width = 13, height = 22)
message('done figures')
