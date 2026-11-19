# Section 2.4 summary figure: gain from each alternative input over native input, by tissue and score, in one panel.
# Estimates (points, 95% CI) = the section 2.4 model in analysis/mcb_analysis.R (input * tissue + (1 | algorithm) +
# (1 | algorithm:input), Kenward-Roger; emmeans contrasts vs native), read from grouping_input_model<suffix>.csv.
# Faint dots behind = the raw gain of each algorithm (same algorithm on the alternative input minus on native input),
# from benchmark_improvement<suffix>.csv; in the 'Mean over tissues' row, each algorithm's mean over the four tissues.
# Nothing is refitted: the estimates are read as saved, and the raw gains are checked against grouping_input_paired<suffix>.csv.
#   Rscript --vanilla analysis/input_summary_figure.R
# Output: figures/mcb/grouping_input_summary<suffix>.{png,pdf}
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(ggplot2)})

# ---------------- config ----------------
suffix     = '_MAP_SAS'
model_set  = 'All matched algorithms'          # as the section 2.4 primary model (6 algorithms incl. BBKNN)
algorithms = c('Harmony', 'BBKNN', 'FastMNN', 'RPCA', 'CCA', 'Scanorama')
highlight  = 'FastMNN'                           # algorithm marked separately among the raw gains
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver',
                  all = 'Mean over tissues')
score_labels  = c(overall = 'Overall', bio = 'Bio-conservation', batch = 'Batch correction')
input_cols    = c(`scVI-latent` = '#7B5EA7', UINMF = '#CCBB44')   # as analysis/benchmark_figures.R / analysis/mcb_analysis.R
dodge         = 0.17                             # vertical offset of the two inputs within a row
out           = paste0('figures/mcb/grouping_input_summary', suffix)
# appendix robustness figure: estimates of the main model and the two refits (sets in grouping_input_model<suffix>.csv)
rob_sets   = c(`All matched algorithms` = 'All six algorithms (main)', `Without BBKNN (embedding only)` = 'Without BBKNN',
               `Without FastMNN` = 'Without fastMNN')
rob_cols   = c('#222222', '#228833', '#EE7733')   # Tol bright; not the species or input colours
rob_shapes = c(21, 22, 24)
rob_dodge  = 0.22
out_rob    = paste0('figures/mcb/grouping_input_robustness', suffix)
# -----------------------------------------

est = read_csv(sprintf('figures/mcb/grouping_input_model%s.csv', suffix), show_col_types = FALSE) |>
  filter(table == 'vs_native', set == model_set) |> select(score, tissue, contrast, est, lo, hi)
stopifnot(nrow(est) == 3 * 5 * 2, setequal(est$tissue, names(tissue_labels)))

imp = read_csv(sprintf('figures/benchmark_summary/benchmark_improvement%s.csv', suffix), show_col_types = FALSE) |>
  mutate(algorithm = sub(' .*', '', method), input = ifelse(input %in% c('scVI-latent', 'UINMF'), input, 'Native')) |>
  filter(algorithm %in% algorithms) |>
  select(algorithm, input, tissue, overall = imp_overall, bio = imp_bio, batch = imp_batch) |>
  pivot_longer(c(overall, bio, batch), names_to = 'score')
stopifnot(nrow(imp) == length(algorithms) * 3 * 4 * 3)
raw = imp |> pivot_wider(names_from = input, values_from = value) |>
  pivot_longer(c(`scVI-latent`, UINMF), names_to = 'contrast', values_to = 'alt') |> mutate(gain = alt - Native)
# check: overall raw gains = grouping_input_paired<suffix>.csv (analysis/mcb_analysis.R section 2.4)
pr = read_csv(sprintf('figures/mcb/grouping_input_paired%s.csv', suffix), show_col_types = FALSE)
chk = raw |> filter(score == 'overall') |> inner_join(pr |> select(algorithm, tissue, contrast, diff), by = c('algorithm', 'tissue', 'contrast'))
stopifnot(nrow(chk) == 48, max(abs(chk$gain - chk$diff)) < 1e-12)
raw = bind_rows(raw, raw |> group_by(algorithm, score, contrast) |> summarise(gain = mean(gain), .groups = 'drop') |> mutate(tissue = 'all'))

rows = rev(names(tissue_labels))                                  # mean at the bottom
pos  = function(t, c) match(t, rows) + ifelse(c == 'scVI-latent', dodge, -dodge)
est = est |> mutate(y = pos(tissue, contrast), score = factor(score, levels = names(score_labels), labels = score_labels))
raw = raw |> mutate(y = pos(tissue, contrast), score = factor(score, levels = names(score_labels), labels = score_labels),
                    hl = ifelse(algorithm == highlight, highlight, 'Other algorithms'))

p = ggplot(est, aes(est, y, colour = contrast)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5, fill = 'grey95') +   # mean row
  geom_vline(xintercept = 0, colour = 'grey30', linewidth = 0.4) +
  geom_point(data = raw, aes(gain, y, shape = hl), size = 1.5, alpha = 0.45, stroke = 0.6) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = 'y', linewidth = 0.9) +
  geom_point(aes(fill = contrast), shape = 21, colour = 'black', size = 2.8, stroke = 0.5) +
  facet_grid(. ~ score, scales = 'free_x') +
  scale_y_continuous(breaks = seq_along(rows), labels = tissue_labels[rows], expand = expansion(add = 0.3)) +
  scale_colour_manual(values = input_cols, name = NULL, labels = \(x) paste(x, '− native')) +
  scale_fill_manual(values = input_cols, guide = 'none') +
  scale_shape_manual(values = setNames(c(17, 16), c(highlight, 'Other algorithms')), name = 'Single algorithm',
                     breaks = c('Other algorithms', highlight)) +
  guides(colour = guide_legend(override.aes = list(shape = 21, fill = input_cols, colour = 'black', size = 2.8, alpha = 1,
                                                   linetype = 1)),
         shape = guide_legend(override.aes = list(colour = 'grey40', alpha = 0.7))) +
  labs(x = 'Gain in improvement over native input (model estimate, 95% CI)', y = NULL) +
  theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.text = element_text(face = 'bold', size = 10.5), strip.background = element_rect(fill = 'grey95'),
        axis.text.y = element_text(colour = 'grey10', face = ifelse(rows == 'all', 'bold', 'plain')),
        legend.position = 'bottom')
ggsave(paste0(out, '.png'), p, width = 10, height = 4.6, dpi = 300, bg = 'white')
ggsave(paste0(out, '.pdf'), p, width = 10, height = 4.6, bg = 'white', device = cairo_pdf)
est |> mutate(across(c(est, lo, hi), \(x) round(x, 3))) |> select(score, tissue, contrast, est, lo, hi) |> as.data.frame() |> print(row.names = FALSE)
cat('wrote', paste0(out, '.{png,pdf}'), '\n')

# ---- appendix: robustness of the estimates to dropping BBKNN or fastMNN ----
rob = read_csv(sprintf('figures/mcb/grouping_input_model%s.csv', suffix), show_col_types = FALSE) |>
  filter(table == 'vs_native', set %in% names(rob_sets)) |> select(set, score, tissue, contrast, est, lo, hi)
stopifnot(nrow(rob) == length(rob_sets) * 3 * 5 * 2)
rob = rob |> mutate(k = match(set, names(rob_sets)), y = match(tissue, rows) + (2 - k) * rob_dodge,
                    set = factor(rob_sets[set], levels = rob_sets),
                    score = factor(score, levels = names(score_labels), labels = score_labels),
                    contrast = factor(paste(contrast, '− native'), levels = paste(names(input_cols), '− native')))
pr2 = ggplot(rob, aes(est, y, colour = set, fill = set, shape = set)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5, fill = 'grey95') +
  geom_vline(xintercept = 0, colour = 'grey30', linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = 'y', linewidth = 0.7) +
  geom_point(size = 2.2, stroke = 0.5, colour = 'black') +
  facet_grid(contrast ~ score, scales = 'free_x') +
  scale_y_continuous(breaks = seq_along(rows), labels = tissue_labels[rows], expand = expansion(add = 0.3)) +
  scale_colour_manual(values = setNames(rob_cols, rob_sets), name = NULL) +
  scale_fill_manual(values = setNames(rob_cols, rob_sets), name = NULL) +
  scale_shape_manual(values = setNames(rob_shapes, rob_sets), name = NULL) +
  labs(x = 'Gain in improvement over native input (model estimate, 95% CI)', y = NULL) +
  theme_bw(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.text = element_text(face = 'bold', size = 10.5), strip.background = element_rect(fill = 'grey95'),
        axis.text.y = element_text(colour = 'grey10', face = ifelse(rows == 'all', 'bold', 'plain')),
        legend.position = 'bottom')
ggsave(paste0(out_rob, '.png'), pr2, width = 10, height = 7.5, dpi = 300, bg = 'white')
ggsave(paste0(out_rob, '.pdf'), pr2, width = 10, height = 7.5, bg = 'white', device = cairo_pdf)
cat('wrote', paste0(out_rob, '.{png,pdf}'), '\n')
