# Bio-conservation score of every integration method under EBCC labels (x) vs under the best-labeller labels (y),
# one panel per tissue and dataset where EBCC is not the best labelling in both species (colon is skipped: EBCC
# best in both). Kendall's tau of the two rankings in each panel. Reads the saved 05 label-sensitivity scores only.
# Usage: Rscript analysis/label_sensitivity_scatter.R  -> figures/cross_tissue/label_sensitivity_bio_scatter.{png,pdf,csv}
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(readr); library(stringr); library(ggplot2) })
setwd(Sys.getenv('PROJECT_ROOT', '.'))

f = \(t) read_csv(sprintf('%s/data/05_evaluated_score/label_sensitivity_scores.csv', t), show_col_types = FALSE)
bm = f('bone_marrow') |> select(dataset, method, bio_ebcc, bio_alt)
long_to_wide = \(d, set_col, alt) d |> select(dataset, method, bio, set = all_of(set_col)) |>
  mutate(set = if_else(set == 'ebcc', 'bio_ebcc', 'bio_alt')) |> pivot_wider(names_from = set, values_from = bio)
d = bind_rows(mutate(bm, tissue = 'bone_marrow'),
              mutate(long_to_wide(f('fetal_liver'), 'label_set'), tissue = 'fetal_liver'),
              mutate(long_to_wide(f('breast_cancer'), 'labels'), tissue = 'breast_cancer')) |>
  filter(!is.na(bio_ebcc), !is.na(bio_alt)) |>
  mutate(dataset = str_remove(dataset, '_mouse$'),
         dataset = factor(dataset, levels = c('full', 'dissected', 'genes2', 'genes4', 'genes6', 'genes8')),
         tissue = factor(tissue, levels = c('bone_marrow', 'fetal_liver', 'breast_cancer')),
         scanvi = if_else(method == 'scANVI', 'scANVI (trained on EBCC labels)', 'Other methods'))
# weighted tau (scipy.stats.weightedtau, hyperbolic weights: swaps among the top methods count more) and its
# one-sided permutation p (9,999 shuffles; smallest possible p = 1e-4), from analysis/tau_extra.py
wt = read_csv('figures/cross_tissue/tau_weighted.csv', show_col_types = FALSE)
taus = wt |> filter(score == 'bio') |> select(tissue, dataset, n = n_methods, tau = weighted_tau, p = weighted_tau_p, tau_plain = tau) |>
  mutate(dataset = factor(dataset, levels = levels(d$dataset)), tissue = factor(tissue, levels = levels(d$tissue)))
# the python scores match these: plain tau recomputed here equals analysis/tau_extra.py's
chk = d |> group_by(tissue, dataset) |> summarise(t = cor(bio_ebcc, bio_alt, method = 'kendall'), .groups = 'drop') |>
  inner_join(taus, by = c('tissue', 'dataset'))
stopifnot(nrow(chk) == 18, max(abs(chk$t - chk$tau_plain)) < 1e-9)
write_csv(left_join(d, taus, by = c('tissue', 'dataset')), 'figures/cross_tissue/label_sensitivity_bio_scatter.csv')
print(as.data.frame(taus))

alt = c(bone_marrow = 'human: Seurat; mouse: EBCC', fetal_liver = 'human: Seurat; mouse: EBCC',
        breast_cancer = 'human: Seurat; mouse: Seurat (cross)')
# permutation p: 1e-4 is the floor (no shuffle reached the observed value), shown as a bound
format_p = \(p) ifelse(p <= 1e-4 + 1e-12, '≤ 1e-4', sprintf('= %.4f', p))
tissue_name = \(t) str_to_sentence(str_replace(t, '_', ' '))
# genes2/4/6/8: the dissected dataset with 20/40/60/80% of shared one-to-one orthologues removed (EBCC_PIPELINE.md 6)
dataset_lab = c(full = 'Full mouse', dissected = 'Dissected mouse',
                genes2 = 'Dissected, 20% orthologues removed', genes4 = 'Dissected, 40% orthologues removed',
                genes6 = 'Dissected, 60% orthologues removed', genes8 = 'Dissected, 80% orthologues removed')
bar_lab = c(full = 'Full', dissected = 'Dissected', genes2 = '−20% orth.', genes4 = '−40% orth.',
            genes6 = '−60% orth.', genes8 = '−80% orth.')
score_lab = c(bio = 'Bio-conservation', overall = 'Overall (0.6 bio + 0.4 batch)')
tau_all = wt |> mutate(score = score_lab[score], tau = weighted_tau)
# sampling-noise range: 5-95% of weighted tau between pairs of bootstrap repeats (dissected dataset, 20,000-cell
# subsamples, EBCC labels; analysis/tau_extra.py), per tissue and score
noise = read_csv('figures/cross_tissue/tau_noise_summary.csv', show_col_types = FALSE) |>
  transmute(tissue, score = score_lab[score], lo = noise_wtau_q05, hi = noise_wtau_q95)
# square panels: each panel's range covers both axes, so y = x is the diagonal
lims = d |> group_by(tissue, dataset) |>
  summarise(lo = min(bio_ebcc, bio_alt), hi = max(bio_ebcc, bio_alt), .groups = 'drop') |>
  mutate(pad = 0.05 * (hi - lo), lo = lo - pad, hi = hi + pad) |>
  pivot_longer(c(lo, hi), values_to = 'v')
theme_row = theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = 'bold', size = 8.5), plot.title = element_text(face = 'bold', size = 12))
row_plot = function(t) {
  dt = filter(d, tissue == t); lt = filter(lims, tissue == t); tt = filter(taus, tissue == t)
  tt = lt |> group_by(dataset) |> summarise(x0 = min(v) + 0.04 * (max(v) - min(v)), .groups = 'drop') |>
    inner_join(tt, by = 'dataset')
  sc = ggplot(dt, aes(bio_ebcc, bio_alt)) +
    geom_blank(data = lt, aes(x = v, y = v)) +
    geom_abline(aes(slope = 1, intercept = 0, linetype = 'y = x (same score under both)'), colour = 'grey55', linewidth = 0.4) +
    geom_point(aes(shape = scanvi, colour = scanvi), size = 1.7, alpha = 0.85) +
    # weighted tau, left-aligned 4% in from the panel's left edge
    geom_text(data = tt, aes(x = x0, y = Inf, label = sprintf('weighted τ = %.2f', tau)),
              hjust = 0, vjust = 1.5, size = 3, fontface = 'bold', inherit.aes = FALSE) +
    facet_wrap(~ dataset, nrow = 1, scales = 'free', labeller = as_labeller(\(x) str_wrap(dataset_lab[x], 20))) +
    scale_shape_manual(values = c('Other methods' = 16, 'scANVI (trained on EBCC labels)' = 17), name = NULL) +
    scale_colour_manual(values = c('Other methods' = 'grey25', 'scANVI (trained on EBCC labels)' = '#C0392B'), name = NULL) +
    scale_linetype_manual(values = 'dashed', name = NULL) +
    scale_x_continuous(n.breaks = 3) + scale_y_continuous(n.breaks = 3) +
    labs(title = sprintf('%s (best labeller: %s)', tissue_name(t), alt[[t]]),
         x = 'Bio-conservation, EBCC labels', y = 'Bio-conservation,\nbest-labeller labels') +
    theme_row + theme(aspect.ratio = 1)
  # dot plot: one row per dataset, bio on the lower lane and overall on the upper lane of each row; the shaded
  # lane behind each dot is that score's sampling-noise range (5-95% between bootstrap subsamples)
  score_cols = c('Bio-conservation' = '#2A7F9E', 'Overall (0.6 bio + 0.4 batch)' = 'grey35')
  bt = filter(tau_all, tissue == t) |>
    mutate(dataset = factor(dataset, levels = rev(names(bar_lab))),
           y = as.numeric(dataset) + if_else(score == score_lab[['bio']], -0.16, 0.16)) |>
    inner_join(filter(noise, tissue == t), by = 'score')
  br = ggplot(bt, aes(tau, y, colour = score)) +
    geom_rect(aes(xmin = lo, xmax = hi, ymin = y - 0.12, ymax = y + 0.12, fill = score,
                  alpha = 'Sampling noise (5–95% between\nbootstrap subsamples)'), colour = NA) +
    geom_point(size = 2.3) +
    geom_text(aes(label = sprintf('%.2f', tau)), hjust = -0.35, size = 2.4, show.legend = FALSE) +
    scale_colour_manual(values = score_cols, name = 'Weighted τ for:') +
    scale_fill_manual(values = score_cols, guide = 'none') +
    scale_alpha_manual(values = 0.18, name = NULL, guide = guide_legend(override.aes = list(fill = 'grey45', colour = NA))) +
    scale_x_continuous(breaks = c(0.6, 0.8, 1), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0.5, 1.02)) +   # clip, don't drop: noise lanes can start below the axis
    scale_y_continuous(breaks = seq_along(levels(bt$dataset)), labels = bar_lab[levels(bt$dataset)],
                       expand = expansion(add = 0.45)) +
    labs(title = ' ', x = 'Weighted Kendall\'s τ', y = NULL) +
    theme_row + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
  patchwork::wrap_plots(sc, br, widths = c(6, 2))
}
p = patchwork::wrap_plots(lapply(levels(d$tissue), row_plot), ncol = 1, guides = 'collect') +
  patchwork::plot_annotation(
    caption = paste('Scatter: one point per integration method; weighted τ: Kendall\'s tau with hyperbolic weights (scipy weightedtau), swaps among the top methods count more.',
                    'Axes differ between panels. Right: weighted \u03c4 for every dataset, bio-conservation (lower dot) and overall score (upper dot).',
                    'Shaded lanes: 5\u201395% of weighted \u03c4 between independent 20,000-cell subsamples of the dissected dataset under EBCC labels (50 bootstrap repeats, full bio score).',
                    'Orth.: shared one-to-one orthologues removed from the dissected dataset.',
                    'Best labeller: per species, the labelling with the highest mean marker AUROC on the selection half (four annotators and EBCC). Colon not shown: EBCC is best in both species.',
                    sep = '\n'),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8.5))) &
  theme(legend.position = 'top')
ggsave('figures/cross_tissue/label_sensitivity_bio_scatter.png', p, width = 16.5, height = 10, dpi = 300)
ggsave('figures/cross_tissue/label_sensitivity_bio_scatter.pdf', p, width = 16.5, height = 10, device = cairo_pdf)
