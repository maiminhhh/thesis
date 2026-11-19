# Does EBCC help where the best single annotator is unsure, and leave it alone where it is sure?
# Run from the repository root. Reads the saved 02 annotation objects; nothing is refitted or rescored.
# Best annotator per dataset and species: highest mean marker AUROC on the SELECTION half among the four
# annotators (figures/cross_tissue/best_annotator_select.csv, from analysis/best_annotator_select.R).
# Evaluation half of cells (as in analysis/vote_vs_ebcc.R), cells the best annotator labels, split into fifths by
# that annotator's own confidence score (Seurat: prediction score labelscore_*; SingleR: margin delta.next,
# labeldelta_*; quintiles within each dataset). Per fifth:
#   1. % of cells where EBCC keeps the best annotator's label
#   2. marker AUROC of EBCC's and of the best annotator's labels on the same cells
#   3. overrides (EBCC label != annotator label, both with markers): % where the cell's marker score
#      ranks higher for EBCC's type than for the annotator's (percentile among evaluation cells, per type)
# plus a histogram of the confidence score. Breast cancer mouse is left out of the figure.
# Usage:  Rscript analysis/seurat_vs_ebcc.R <tissue> <species>   -> figures/cross_tissue/seurat_vs_ebcc/<tissue>_<species>.csv
#         Rscript analysis/seurat_vs_ebcc.R                      figure

# helpers from analysis/vote_vs_ebcc.R (and through it analysis/visualisations.R)
# (Seurat-free, as analysis/ebcc_support.R: Seurat no longer loads under the r45 env; the figure step reads CSVs only)
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
vis = parse('analysis/visualisations.R'); vis_h = sapply(vis, \(e) deparse(e)[1])
for (i in seq_len(which(startsWith(vis_h, 'for (tissue in names(tissues))')) - 1))
  if (!grepl('^(source|library)\\(|run_only', vis_h[i])) eval(vis[[i]], globalenv())
ex = parse('analysis/vote_vs_ebcc.R'); ex_h = sapply(ex, \(e) deparse(e)[1])
for (i in (which(startsWith(ex_h, 'for (e in ex'))[1] + 1):(which(startsWith(ex_h, 'args = commandArgs')) - 1))
  eval(ex[[i]], globalenv())

out_dir = 'figures/cross_tissue/seurat_vs_ebcc'
ann_keys = c('Seurat' = 'seurat', 'SingleR' = 'singleR', 'Seurat (cross)' = 'seurat_cross', 'SingleR (cross)' = 'singleR_cross')
score_cols = c(seurat = 'labelscore_seurat', singleR = 'labeldelta_singleR',
               seurat_cross = 'labelscore_seurat_cross', singleR_cross = 'labeldelta_singleR_cross')
excluded = tibble(tissue = 'breast_cancer', species = 'mouse')
best_ann = read_csv('figures/cross_tissue/best_annotator_select.csv', show_col_types = FALSE) |>
  filter(best) |> select(tissue, species, annotator)

compute_one = function(tissue, species) {
  t0 = Sys.time()
  tick = function(msg) message(sprintf('[%s %s %5.1f min] %s', tissue, species,
                                       as.numeric(Sys.time() - t0, units = 'mins'), msg))
  annotator = best_ann$annotator[best_ann$tissue == tissue & best_ann$species == species]
  k = ann_keys[[annotator]]
  is_seurat = startsWith(k, 'seurat')
  a = tissues[[tissue]]$annotation[[species]]
  obj = readRDS(a$obj)
  assay = DefaultAssay(obj)
  rownames(obj) = obj[[assay]]@meta.data$gene_name
  X = LayerData(obj, assay = assay, layer = 'data')
  meta = obj@meta.data
  rm(obj); gc(); tick(sprintf('object loaded; best annotator %s', annotator))
  set.seed(1)
  held = !seq_len(ncol(X)) %in% sample(ncol(X), floor(ncol(X) / 2))
  markers = lapply(annotation_markers(a), intersect, rownames(X))
  markers = markers[lengths(markers) > 0]
  score = sapply(markers, \(g) Matrix::colMeans(X[g, held, drop = FALSE]))   # held cells x types
  colnames(score) = key(colnames(score))
  rm(X); gc()
  pct = apply(score, 2, \(x) rank(x) / length(x))                            # per-type percentile
  ch = read_csv(file.path(dirname(a$obj), str_glue('{species}_ebcc_frac_auroc_chosen.csv')), show_col_types = FALSE)

  ref = key(as.character(meta[[annotators[[k]]]])[held])
  ebcc = key(as.character(meta$ebcc_label)[held])
  s_score = meta[[score_cols[[k]]]][held]
  keep = !is.na(ref) & !is.na(s_score)                                       # SingleR: pruned calls left out
  bin = rep(NA_integer_, length(ref))
  bin[keep] = dplyr::ntile(s_score[keep], 5)
  mean_auc_in = function(l, idx) {
    types = colnames(score)[vapply(colnames(score), \(k) { p = l[idx] == k; any(p) && !all(p) }, logical(1))]
    mean(vapply(types, \(k) auroc(score[idx, k], l[idx] == k), numeric(1)))
  }
  res = purrr::map_dfr(1:5, function(b) {
    idx = which(bin == b)
    ov = idx[ref[idx] != ebcc[idx]]
    scored = ov[ref[ov] %in% colnames(score) & ebcc[ov] %in% colnames(score)]
    pe = pct[cbind(scored, match(ebcc[scored], colnames(score)))]; pr = pct[cbind(scored, match(ref[scored], colnames(score)))]
    win = pe > pr
    ci = if (length(win)) prop.test(sum(win), length(win))$conf.int else c(NA, NA)
    tibble(bin = b, n = length(idx), score_min = min(s_score[idx]), score_max = max(s_score[idx]),
           agree = mean(ref[idx] == ebcc[idx]), n_override = length(ov), n_override_scored = length(scored),
           override_win = if (length(win)) mean(win) else NA_real_, win_lo = ci[1], win_hi = ci[2],
           auroc_ebcc = mean_auc_in(ebcc, idx), auroc_best = mean_auc_in(ref, idx),
           # every override, for the net-effect panel: EBCC's type fits better / worse / equally (tied percentile,
           # e.g. both marker scores 0), or one of the two types has no markers in the panel (unscored)
           n_fix = sum(pe > pr), n_break = sum(pe < pr), n_tie = sum(pe == pr), n_unscored = length(ov) - length(scored))
  }) |> mutate(tissue = tissue, species = species, annotator = annotator,
               score_name = if (is_seurat) 'prediction score' else 'delta', chosen_frac = ch$frac,
               n_held = sum(held), n_used = sum(keep), .before = 1)
  # histogram of the score on the same cells: Seurat bin width 0.05 on [0, 1] (a score of 1 in the last bin);
  # SingleR delta 20 equal bins on [0, max]
  s = s_score[keep]
  width = if (is_seurat) 0.05 else max(s) / 20
  hist = tibble(bin_left = width * pmin(floor(s / width + 1e-9), 19)) |> dplyr::count(bin_left) |>
    mutate(tissue = tissue, species = species, annotator = annotator, bin_width = width, pct = n / length(s),
           pct_one = if (is_seurat) mean(s >= 1 - 1e-9) else NA_real_, .before = 1)
  tick('done')
  list(res = res, hist = hist)
}

args = commandArgs(TRUE)
if (length(args) == 2) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  r = compute_one(args[1], args[2])
  write_csv(r$res, file.path(out_dir, sprintf('%s_%s.csv', args[1], args[2])))
  write_csv(r$hist, file.path(out_dir, sprintf('%s_%s_hist.csv', args[1], args[2])))
  quit(save = 'no')
}

# ---- figure: one column per tissue; A: score distribution (human row, mouse row), B-C: by quintile ----
# (a marker-AUROC-by-quintile panel, p3, is built but not shown: it weights cell types equally, so a few cells moved
#  into a rare type move it a lot; C counts cells instead)
datasets = expand_grid(tissue = tissue_order, species = c('human', 'mouse')) |> anti_join(excluded, by = c('tissue', 'species'))
read_part = \(suffix) purrr::pmap_dfr(datasets, \(tissue, species)
  read_csv(file.path(out_dir, sprintf('%s_%s%s.csv', tissue, species, suffix)), show_col_types = FALSE))
res = read_part('') |> mutate(tissue = factor(tissue, levels = tissue_order))
hist = read_part('_hist') |> mutate(tissue = factor(tissue, levels = tissue_order))
cuts = res |> filter(bin < 5) |> select(tissue, species, annotator, bin, cut = score_max)
quintiles = cuts |> mutate(bin = paste0('cut_', bin)) |> pivot_wider(names_from = bin, values_from = cut) |>
  left_join(distinct(res, tissue, species, score_name, n_held, n_used), by = c('tissue', 'species')) |>
  left_join(distinct(hist, tissue, species, pct_one), by = c('tissue', 'species'))
write_csv(res, 'figures/cross_tissue/seurat_vs_ebcc.csv')
write_csv(hist, 'figures/cross_tissue/seurat_vs_ebcc_score_hist.csv')
write_csv(quintiles, 'figures/cross_tissue/seurat_vs_ebcc_quintiles.csv')

species_vals = c(human = species_cols[['hsapiens']], mouse = species_cols[['mmusculus']])
tissue_name = \(t) str_to_sentence(str_replace(t, '_', ' '))
# panel A facets: human row then mouse row, one column per tissue; the excluded dataset keeps an empty slot
facet_levels = expand_grid(species = c('human', 'mouse'), tissue = tissue_order) |>
  left_join(distinct(res, tissue = as.character(tissue), species, annotator, score_name), by = c('tissue', 'species')) |>
  mutate(id = paste(tissue, species),
         lab = if_else(is.na(annotator), sprintf('%s, %s\nexcluded', tissue_name(tissue), species),
                       sprintf('%s, %s\n%s (%s)', tissue_name(tissue), species, annotator,
                               if_else(score_name == 'delta', '\u0394', 'prediction score'))))
facet_a = \(d) mutate(d, facet = factor(paste(tissue, species), levels = facet_levels$id))
ones = distinct(hist, tissue, species, pct_one) |> filter(!is.na(pct_one)) |>
  mutate(label = sprintf('%.0f%% at 1', 100 * pct_one)) |> facet_a()
p1 = ggplot(facet_a(hist)) +
  geom_rect(aes(xmin = bin_left, xmax = bin_left + bin_width, ymin = 0, ymax = pct, fill = species), alpha = 0.6) +
  geom_vline(data = facet_a(cuts), aes(xintercept = cut, colour = species), linetype = 'dashed', linewidth = 0.4) +
  geom_text(data = ones, aes(x = -Inf, y = Inf, label = label, colour = species), hjust = -0.1, vjust = 1.5, size = 3) +
  facet_wrap(~ facet, nrow = 2, scales = 'free', drop = FALSE, labeller = as_labeller(setNames(facet_levels$lab, facet_levels$id))) +
  scale_fill_manual(values = species_vals, guide = 'none') + scale_colour_manual(values = species_vals, guide = 'none') +
  scale_y_continuous(labels = scales::percent) +
  labs(x = 'Best annotator\'s confidence score', y = '% of cells') +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = 'bold', hjust = 0, size = 9))
dodge = position_dodge(width = 0.35)
row_panel = function(p, ylab, top = FALSE) p +
  facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(tissue_name)) +
  scale_colour_manual(values = species_vals, name = NULL) +
  scale_x_continuous(breaks = 1:5) +
  labs(x = NULL, y = ylab) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = if (top) element_text(face = 'bold', hjust = 0) else element_blank())
p2 = row_panel(ggplot(res, aes(bin, agree, colour = species)) + geom_line() + geom_point(size = 1.8) +
                 scale_y_continuous(labels = scales::percent, limits = c(0, 1)),
               'EBCC keeps the best\nannotator\'s label', top = TRUE) + guides(colour = 'none')
p3 = row_panel(ggplot(res, aes(bin, colour = species)) +
                 geom_line(aes(y = auroc_ebcc, linetype = 'EBCC (chosen f)')) + geom_point(aes(y = auroc_ebcc), size = 1.8) +
                 geom_line(aes(y = auroc_best, linetype = 'Best annotator')) +
                 scale_linetype_manual(values = c('EBCC (chosen f)' = 'solid', 'Best annotator' = 'dashed'), name = 'C: labels scored'),
               'Mean marker\nAUROC') +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2))
# C: EBCC's overrides per 100 cells of the fifth: supported by cell markers (EBCC's type fits the markers better) up,
# not supported (the annotator's type fits better) down, unscored overrides (a type without markers) stacked on the
# supported ones in grey; line = net effect (n_fix / n_break below = supported / not supported)
eff = res |> mutate(fix = 100 * n_fix / n, brk = -100 * n_break / n, uns = 100 * n_unscored / n, net = fix + brk)
dodge_bar = position_dodge(width = 0.8)
p4 = row_panel(ggplot(eff, aes(bin, colour = species)) +
                 geom_hline(yintercept = 0, colour = 'grey40', linewidth = 0.4) +
                 geom_col(aes(y = fix + uns, group = species, fill = 'Override not scorable (type without markers)'),
                          position = dodge_bar, width = 0.75, colour = NA) +
                 geom_col(aes(y = fix, group = species, fill = species), position = dodge_bar, width = 0.75, colour = NA, alpha = 0.45) +
                 geom_col(aes(y = brk, group = species, fill = species), position = dodge_bar, width = 0.75, colour = NA, alpha = 0.45) +
                 geom_line(aes(y = net, group = species), position = position_dodge(width = 0.8), linewidth = 0.7) +
                 geom_point(aes(y = net, group = species, shape = 'Net effect (supported − not supported)'), position = position_dodge(width = 0.8), size = 1.9) +
                 scale_fill_manual(values = c(species_vals, 'Override not scorable (type without markers)' = 'grey80'),
                                   breaks = 'Override not scorable (type without markers)', name = NULL) +
                 scale_shape_manual(values = 16, name = 'C:') +
                 scale_y_continuous(labels = \(x) sprintf('%+g', x)),
               'EBCC\'s changes per 100 cells:\nsupported by markers (+)\nnot supported (−)') +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16)), shape = guide_legend(order = 2)) +
  labs(x = 'Best annotator\'s confidence score, quintile within dataset: 1 = least confident 20% of cells, 5 = most confident')
p = patchwork::wrap_plots(p1, p2, p4, ncol = 1, guides = 'collect', heights = c(2, 1, 1.15)) +
  patchwork::plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'top')
writeLines(paste('Evaluation half of cells (not used to choose f); breast cancer mouse excluded. Best annotator: highest mean marker AUROC on the selection half',
                    'among the four annotators (named in A); EBCC uses its labels, not its scores. Cells it leaves unlabelled (SingleR pruned calls) are excluded.',
                    'A: distribution of the best annotator\'s confidence score (Seurat: prediction score; SingleR: \u0394, margin to the next-best label);',
                    'dashed lines mark the quintile boundaries used in B–C. Many Seurat scores are exactly 1, so the top quintiles split tied cells.',
                    'C: cells where EBCC changes the annotator\'s label, per 100 cells of the fifth: supported by cell markers if the cell\'s marker score is at a higher percentile',
                    '(among evaluation cells) for EBCC\'s type than for the annotator\'s, not supported if lower (no ties occurred); line: net effect (supported minus not supported).',
                    sep = '\n'), 'figures/cross_tissue/seurat_vs_ebcc_caption.txt')   # caption for the thesis text, not drawn
save_plot(p, 'figures/cross_tissue/seurat_vs_ebcc.png', width = 13, height = 11)
ggsave('figures/cross_tissue/seurat_vs_ebcc.pdf', p, width = 13, height = 11, device = cairo_pdf)
message('done figure')
