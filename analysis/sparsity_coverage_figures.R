# Two figures for the annotation-disagreement discussion (reads saved CSVs only).
# 1. sparsity_concentration.png (analysis/check_sparsity_concentration.R): within each label of one annotator, cells in fifths of
#    panel genes detected; share the OTHER annotator puts into its own dominant label, for SingleR and for Seurat.
# 2. singler_gene_detection.png (analysis/check_singler_genes.R): per reference cell type, how often the genes SingleR uses to
#    recognise it are detected in the query (x) against how many more cells SingleR gives the type than Seurat (y).
# Breast cancer mouse left blank in both. Usage: Rscript45 analysis/sparsity_coverage_figures.R
setwd(Sys.getenv('PROJECT_ROOT', '.')); source('library.R')
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2) })
tissue_order = c('colon', 'bone_marrow', 'fetal_liver', 'breast_cancer')
tissue_name = \(t) str_to_sentence(str_replace(t, '_', ' '))
panel_levels = as.vector(t(outer(tissue_name(tissue_order), c('human', 'mouse'), \(a, b) sprintf('%s, %s', a, b))))
panel_of = \(t, s) factor(sprintf('%s, %s', tissue_name(t), s), levels = panel_levels)
ann_cols = c(SingleR = '#D2691E', Seurat = '#2A7F9E')
theme_fig = theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = 'bold', hjust = 0), legend.position = 'top')

# ---- 1. sparsity vs concentration ----
sc = read_csv('figures/cross_tissue/sparsity_concentration.csv', show_col_types = FALSE) |>
  mutate(annotator = if_else(str_starts(direction, 'SingleR'), 'SingleR', 'Seurat'),
         panel = panel_of(tissue, species))
doms = sc |> distinct(panel, annotator, dominant) |> pivot_wider(names_from = annotator, values_from = dominant) |>
  mutate(label = sprintf('SingleR → %s\nSeurat → %s', SingleR, Seurat))
p1 = ggplot(sc, aes(q, to_dominant, colour = annotator)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_text(data = doms, aes(x = 5.15, y = Inf, label = label), inherit.aes = FALSE, hjust = 1, vjust = 1.2,
            size = 2.6, lineheight = 0.9, colour = 'grey25') +
  facet_wrap(~ panel, ncol = 2, drop = FALSE) +
  scale_colour_manual(values = ann_cols, name = 'Annotator assigning its dominant label:') +
  scale_x_continuous(breaks = 1:5, labels = c('1\nfewest', '2', '3', '4', '5\nmost'), limits = c(0.85, 5.15)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = 'Fifth of panel genes detected per cell (within the other annotator\'s label)',
       y = 'Cells given the annotator\'s dominant label') +
  theme_fig
ggsave('figures/cross_tissue/sparsity_concentration.png', p1, width = 9, height = 10, dpi = 300)
ggsave('figures/cross_tissue/sparsity_concentration.pdf', p1, width = 9, height = 10, device = cairo_pdf)
writeLines(paste('Same-species reference. For each annotator, the cells the other annotator does not assign to this annotator\'s dominant label',
                 '(its most frequent label, named in each panel) are split into fifths by the number of panel genes detected (count > 0),',
                 'within each of the other annotator\'s labels with at least 200 cells, so cell type is roughly held fixed.',
                 'y: share of those cells the annotator assigns to its dominant label. SingleR pruned calls excluded.',
                 'Breast cancer mouse left blank: its reference has no Cancer Epithelial type.', sep = '\n'),
           'figures/cross_tissue/sparsity_concentration_caption.txt')

# ---- 2. detection of SingleR's genes vs over-assignment ----
sg = list.files('figures/cross_tissue/singler_genes', '\\.csv$', full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |> bind_rows() |>
  mutate(panel = panel_of(tissue, species), diff = share_singler - share_seurat) |>
  group_by(panel) |>
  mutate(top_detect = mean_detect == max(mean_detect), dominant = share_singler == max(share_singler),
         label = if_else(abs(diff) >= 0.05 | dominant | top_detect, type, NA_character_)) |> ungroup()
p2 = ggplot(sg, aes(mean_detect, diff)) +
  geom_hline(data = distinct(sg, panel), aes(yintercept = 0), colour = 'grey60', linewidth = 0.4) +
  geom_point(aes(size = share_singler, fill = dominant), shape = 21, colour = 'grey20', stroke = 0.4, alpha = 0.9) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.5, min.segment.length = 0.2, seed = 1, max.overlaps = 20,
                           colour = 'grey20', box.padding = 0.35) +
  facet_wrap(~ panel, ncol = 2, drop = FALSE, scales = 'free_x') +
  scale_fill_manual(values = c(`TRUE` = '#D2691E', `FALSE` = 'grey80'), labels = c(`TRUE` = 'SingleR\'s dominant label', `FALSE` = 'Other types'),
                    name = NULL) +
  scale_size_area(max_size = 6, labels = scales::percent, breaks = c(0.05, 0.25, 0.5), name = 'Cells SingleR assigns') +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = \(x) sprintf('%+.0f', 100 * x)) +
  labs(x = 'Detection of the type\'s SingleR genes in the query (mean share of cells detecting each gene)',
       y = 'Cells assigned by SingleR minus by Seurat (percentage points)') +
  theme_fig + guides(fill = guide_legend(override.aes = list(size = 4)))
ggsave('figures/cross_tissue/singler_gene_detection.png', p2, width = 10, height = 11, dpi = 300)
ggsave('figures/cross_tissue/singler_gene_detection.pdf', p2, width = 10, height = 11, device = cairo_pdf)
writeLines(paste('One point per reference cell type (same-species reference). SingleR genes: the genes SingleR selects as up-regulated in the type',
                 'in any pairwise comparison with another type (retrained exactly as in the annotation: de.method = t, genes shared with the',
                 'query; the retrained model reproduces the saved SingleR labels for 100% of 3,000 sampled cells per dataset).',
                 'x: mean over those genes of the share of query cells detecting the gene (count > 0). y: share of cells SingleR assigns the type',
                 'minus the share Seurat assigns it (cells both annotators label). Labelled: types with a difference of at least 5 points, the',
                 'dominant SingleR label and the type with the most detected genes. Breast cancer mouse left blank: its reference has no Cancer',
                 'Epithelial type.', sep = '\n'),
           'figures/cross_tissue/singler_gene_detection_caption.txt')
cat('done figures\n')
