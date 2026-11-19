# Is the chosen annotator really the best? Per dataset and species, the four annotators' mean marker AUROC on
# the SELECTION half (the values that pick the best annotator; analysis/best_annotator_select.R) and on the held-out
# EVALUATION half with 95% bootstrap intervals (analysis/ebcc_support.R, 'all' cells; SingleR on the cells it labels).
# Reads saved CSVs only. Usage: Rscript analysis/best_annotator_plot.R  -> figures/cross_tissue/best_annotator.{png,pdf,csv}
# helpers without Seurat (it no longer loads under the r45 env): analysis/visualisations.R config/helpers minus its
# source()/library() calls, then analysis/vote_vs_ebcc.R's helpers (as analysis/ebcc_support.R)
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
vis = parse('analysis/visualisations.R'); vis_h = sapply(vis, \(e) deparse(e)[1])
for (i in seq_len(which(startsWith(vis_h, 'for (tissue in names(tissues))')) - 1))
  if (!grepl('^(source|library)\\(|run_only', vis_h[i])) eval(vis[[i]], globalenv())
ex = parse('analysis/vote_vs_ebcc.R'); ex_h = sapply(ex, \(e) deparse(e)[1])
for (i in (which(startsWith(ex_h, 'for (e in ex'))[1] + 1):(which(startsWith(ex_h, 'args = commandArgs')) - 1))
  eval(ex[[i]], globalenv())

ann_levels = c('Seurat', 'SingleR', 'Seurat (cross)', 'SingleR (cross)')
ebcc = 'EBCC (chosen f)'
all_levels = c(ann_levels, ebcc)
sel = read_csv('figures/cross_tissue/best_annotator_select.csv', show_col_types = FALSE) |>
  select(tissue, species, annotator, auroc_select = mean_auroc, best)
# EBCC for comparison (not a candidate): selection-half AUROC at the chosen f, from the frac grid that chose f
grid_col = \(g) intersect(c('mean_auroc_selection', 'mean_auroc_select'), names(g))[1]
sel_ebcc = map_dfr(tissue_order, \(t) map_dfr(c('human', 'mouse'), function(sp) {
  g = read_csv(sprintf('%s/data/02_annotated_data/%s_ebcc_frac_auroc.csv', t, sp), show_col_types = FALSE)
  tibble(tissue = t, species = sp, annotator = ebcc, auroc_select = g[[grid_col(g)]][as.logical(g$chosen)], best = FALSE)
}))
ev = read_csv('figures/cross_tissue/ebcc_hard_cells.csv', show_col_types = FALSE) |>
  filter(set == 'all', labelling %in% all_levels) |>
  select(tissue, species, annotator = labelling, auroc_eval = mean_auroc, lo, hi)
d = inner_join(bind_rows(sel, sel_ebcc), ev, by = c('tissue', 'species', 'annotator')) |>
  group_by(tissue, species) |>
  # best labelling of all five (four annotators and EBCC), chosen on the selection half
  mutate(best_annotator = best, best = auroc_select == max(auroc_select),
         rank_select = rank(-auroc_select), rank_eval = rank(-auroc_eval),
         lead_select = auroc_select[best] - max(auroc_select[!best]),               # best vs runner-up of the five
         lead_eval = auroc_eval[best] - max(auroc_eval[!best])) |>
  ungroup() |>
  mutate(tissue = factor(tissue, levels = tissue_order), annotator = factor(annotator, levels = rev(all_levels)),
         panel = factor(sprintf('%s, %s', str_to_sentence(str_replace(tissue, '_', ' ')), species),
                        levels = as.vector(t(outer(str_to_sentence(str_replace(tissue_order, '_', ' ')),
                                                   c('human', 'mouse'), \(a, b) sprintf('%s, %s', a, b))))))
write_csv(d, 'figures/cross_tissue/best_annotator.csv')

leads = d |> filter(best) |>
  mutate(label = sprintf('lead over runner-up\nselection %+.3f\nevaluation %+.3f', lead_select, lead_eval))
half_cols = c('Selection half (used to choose f)' = 'grey15', 'Evaluation half (held out, 95% CI)' = '#2A7F9E')
p = ggplot(d, aes(y = annotator)) +
  # row shading = rank on the held-out evaluation half (1 = highest, darkest), as the reported comparison
  geom_rect(aes(ymin = as.numeric(annotator) - 0.45, ymax = as.numeric(annotator) + 0.45, fill = factor(rank_eval)),
            xmin = -Inf, xmax = Inf) +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = 'Evaluation half (held out, 95% CI)'),
                 height = 0, linewidth = 0.8, position = position_nudge(y = -0.15)) +
  geom_point(aes(x = auroc_eval, colour = 'Evaluation half (held out, 95% CI)'), size = 2.4,
             position = position_nudge(y = -0.15)) +
  geom_point(aes(x = auroc_select, colour = 'Selection half (used to choose f)'), shape = 21, fill = 'white',
             size = 2.4, stroke = 1, position = position_nudge(y = 0.15)) +
  facet_wrap(~ panel, nrow = 2, dir = 'v') +
  scale_colour_manual(values = half_cols, name = NULL) +
  # warm (opposite the teal points), even lightness steps; every step keeps >= 3:1 contrast with the teal points
  scale_fill_manual(values = c('1' = '#F6CDB0', '2' = '#F8DAC4', '3' = '#FAE5D6', '4' = '#FCF0E6', '5' = '#FEF8F5'),
                    name = 'Rank on evaluation half', guide = guide_legend(nrow = 1, order = 2)) +
  scale_x_continuous(breaks = seq(0.5, 1, 0.1)) +
  scale_y_discrete(limits = rev(all_levels), expand = expansion(add = 0.5)) +
  labs(x = 'Mean marker AUROC (mean over the cell types each annotator calls)', y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'top', panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = 'bold', hjust = 0),
        plot.caption = element_text(hjust = 0, size = 9))
# caption for the thesis text, not drawn
writeLines(paste('Row shading: rank of the four annotators and EBCC (chosen f) by mean marker AUROC on the held-out evaluation half; 1 (darkest) = highest.',
                 'Open circles: the selection half (the half used to choose f), for reference.',
                 'Filled points: the held-out evaluation half; bars: 95% bootstrap intervals (200 resamples of cells).',
                 'SingleR scored on the cells it labels (pruned calls excluded).',
                 'Lead of the chosen labelling over the runner-up on each half: best_annotator.csv (lead_select, lead_eval).', sep = '\n'),
           'figures/cross_tissue/best_annotator_caption.txt')
save_plot(p, 'figures/cross_tissue/best_annotator.png', width = 13, height = 8)
ggsave('figures/cross_tissue/best_annotator.pdf', p, width = 13, height = 8)
print(as.data.frame(d |> filter(best) |> select(tissue, species, annotator, auroc_select, auroc_eval, rank_eval, lead_select, lead_eval)))

# ---- bar version, for comparison: evaluation half as bars from chance (0.5), selection half as a tick ----
pb = ggplot(d, aes(y = annotator)) +
  geom_rect(aes(xmin = 0.5, xmax = auroc_eval, ymin = as.numeric(annotator) - 0.35, ymax = as.numeric(annotator) + 0.35,
                fill = ifelse(best, 'Chosen annotator', 'Other annotators'))) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.5, colour = 'grey15') +
  geom_point(aes(x = auroc_select, shape = 'Selection half (used to choose f)'), size = 6, colour = '#E8A33D') +
  facet_wrap(~ panel, nrow = 2, dir = 'v') +
  scale_fill_manual(values = c('Chosen annotator' = '#2A7F9E', 'Other annotators' = 'grey75'), name = 'Evaluation half (held out):') +
  scale_shape_manual(values = c('Selection half (used to choose f)' = 124), name = NULL) +
  scale_y_discrete(limits = rev(ann_levels)) +
  scale_x_continuous(limits = c(0.5, 1), breaks = seq(0.5, 0.9, 0.1), expand = c(0, 0)) +
  labs(x = 'Mean marker AUROC (bars start at 0.5, chance level)', y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'top', panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = 'bold', hjust = 0),
        plot.caption = element_text(hjust = 0, size = 9), panel.spacing.x = unit(1.4, 'lines'))
writeLines(paste('Bars: held-out evaluation half, with 95% bootstrap intervals (200 resamples of cells). Tick: selection half, the values that choose the annotator.',
                 'SingleR scored on the cells it labels (pruned calls excluded).', sep = '\n'),
           'figures/cross_tissue/best_annotator_bars_caption.txt')
save_plot(pb, 'figures/cross_tissue/best_annotator_bars.png', width = 13, height = 7)
ggsave('figures/cross_tissue/best_annotator_bars.pdf', pb, width = 13, height = 7)
