# Appendix: is the best single annotator's own confidence score consistent with marker AUROC?
# Reads figures/cross_tissue/seurat_vs_ebcc.csv (analysis/seurat_vs_ebcc.R; nothing recomputed): evaluation half, cells the best
# annotator labels, split into quintiles of that annotator's confidence score within each dataset (Seurat: prediction score;
# SingleR: delta, margin to the next-best label); auroc_best = mean marker AUROC of the best annotator's labels on the cells
# of the quintile (types present in and not covering the whole quintile, each weighted equally). Breast cancer mouse excluded.
#   Rscript --vanilla analysis/appendix_confidence_auroc.R
# Output: figures/cross_tissue/confidence_auroc.{png,pdf}, confidence_auroc_caption.txt; printout (values, Spearman per dataset).
suppressPackageStartupMessages({library(dplyr); library(readr); library(tidyr); library(ggplot2); library(ggrepel)})

tissue_labels  = c(colon = 'Colon', bone_marrow = 'Bone marrow', fetal_liver = 'Fetal liver', breast_cancer = 'Breast cancer')
species_cols   = c(human = '#4477AA', mouse = '#EE6677')   # as analysis/visualisations.R
out = 'figures/cross_tissue/confidence_auroc'

r = read_csv('figures/cross_tissue/seurat_vs_ebcc.csv', show_col_types = FALSE) |>
  mutate(tissue = factor(tissue, levels = names(tissue_labels)), species = factor(species, levels = c('mouse', 'human')),
         score = ifelse(score_name == 'delta', 'Δ', 'prediction score'), lab = sprintf('%s, %s', annotator, score))
stopifnot(nrow(r) == 7 * 5)
rho = r |> group_by(tissue, species, lab) |> summarise(rho = cor(bin, auroc_best, method = 'spearman'),
                                                       q1 = auroc_best[bin == 1], q5 = auroc_best[bin == 5], .groups = 'drop')
# annotator named in the bottom-right corner of each panel, mouse above human, in the species colour
ends = r |> filter(bin == 5) |> mutate(label = sprintf('%s: %s (%s)', ifelse(species == 'mouse', 'Mouse', 'Human'), annotator, score),
                                       y = ifelse(species == 'mouse' & tissue != 'breast_cancer', 0.625, 0.595))

p = ggplot(r, aes(bin, auroc_best, colour = species)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_text(data = ends, aes(x = 5.15, y = y, label = label), size = 2.9, hjust = 1, show.legend = FALSE) +
  facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(tissue_labels)) +
  scale_colour_manual(values = species_cols, labels = c(human = 'Human', mouse = 'Mouse'), name = NULL, breaks = c('mouse', 'human')) +
  scale_x_continuous(breaks = 1:5, limits = c(0.85, 5.15), expand = expansion(add = 0.1)) +
  scale_y_continuous(limits = c(0.58, 1), breaks = seq(0.6, 1, 0.1)) +
  labs(x = 'Confidence quintile of the best annotator (1 = least confident 20% of cells, 5 = most confident)',
       y = 'Mean marker AUROC of\nthe best annotator\'s labels') +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = 'bold', size = 10),
        legend.position = 'top')
ggsave(paste0(out, '.png'), p, width = 10, height = 3.8, dpi = 300, bg = 'white')
ggsave(paste0(out, '.pdf'), p, width = 10, height = 3.8, device = cairo_pdf, bg = 'white')

writeLines(paste(
  'Marker AUROC of the best single annotator\'s labels by that annotator\'s own confidence. Evaluation half of cells (not used to choose f);',
  'cells the best annotator labels, split into fifths of its confidence score within each dataset (Seurat label transfer: prediction score;',
  'SingleR: Δ, the margin to the next-best label). Points: mean marker AUROC of the annotator\'s labels on the cells of each fifth, averaged',
  'over the cell types it assigns there. The best annotator (named at the end of each line) is the one with the highest mean marker AUROC on',
  'the selection half; quintiles as in Figure X. Breast cancer mouse excluded.'), paste0(out, '_caption.txt'))

options(width = 200)
r |> select(tissue, species, lab, bin, auroc_best) |> mutate(auroc_best = round(auroc_best, 3)) |>
  pivot_wider(names_from = bin, values_from = auroc_best, names_prefix = 'q') |> as.data.frame() |> print(row.names = FALSE)
cat('\nSpearman rho (quintile vs AUROC) per dataset:\n'); print(as.data.frame(rho |> mutate(across(where(is.double), \(x) round(x, 3)))), row.names = FALSE)
cat(sprintf('Q1 -> Q5 gain: median %.3f, range %.3f to %.3f\n', median(rho$q5 - rho$q1), min(rho$q5 - rho$q1), max(rho$q5 - rho$q1)))
