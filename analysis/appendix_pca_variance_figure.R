# Appendix: cumulative explained variance (top) and scree (bottom) of each species' own PCA, per tissue,
# dissected mouse case, evaluation cells. Data from analysis/appendix_pca_variance_compute.py (run that first).
# No footnotes in the image; caption in pca_variance_caption.txt.
# Usage: R_LIBS_SITE= R_LIBS_USER= <r45>/bin/Rscript --vanilla analysis/appendix_pca_variance_figure.R
#   -> figures/benchmark_summary/pca_variance.{png,pdf}, pca_variance_summary.csv, pca_variance_caption.txt
setwd(Sys.getenv('PROJECT_ROOT', '.'))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

# ---- config ----
in_file = 'figures/benchmark_summary/pca_variance.csv'
out = 'figures/benchmark_summary/pca_variance'
tissue_labels = c(colon = 'Colon', bone_marrow = 'Bone marrow', fetal_liver = 'Fetal liver', breast_cancer = 'Breast cancer')
species_cols = c(human = '#4477AA', mouse = '#EE6677')        # as the UMAP figures
ref_pcs = c(10, 20, 50)                                       # dotted guides: PC counts used in the pipeline
scree_pcs = 30                                                # scree: first 30 PCs

d = read.csv(in_file) |>
  mutate(species = recode(species, hsapiens = 'human', mmusculus = 'mouse'),
         tissue_lab = factor(sprintf('%s (%d genes)', tissue_labels[tissue], n_genes),
                             levels = unique(sprintf('%s (%d genes)', tissue_labels[tissue], n_genes))[
                               order(match(unique(tissue), names(tissue_labels)))]))
summ = d |> group_by(tissue, species, n_cells, n_genes) |>
  summarise(PC1 = 100 * var_ratio[1], first10 = 100 * sum(var_ratio[pc <= 10]), first20 = 100 * sum(var_ratio[pc <= 20]),
            first50 = 100 * sum(var_ratio[pc <= 50]), pcs_for_80 = min(pc[cum_ratio >= 0.8]),
            pcs_for_90 = min(pc[cum_ratio >= 0.9]), .groups = 'drop') |>
  arrange(match(tissue, names(tissue_labels)), species)
write.csv(summ, paste0(out, '_summary.csv'), row.names = FALSE)

theme_pc = theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.spacing.x = unit(1.2, 'lines'), strip.background = element_blank(),
        strip.text = element_text(face = 'bold', size = 12), legend.position = 'bottom')
guides_df = expand.grid(tissue_lab = levels(d$tissue_lab), x = ref_pcs) |>
  left_join(distinct(d, tissue_lab, n_genes), by = 'tissue_lab') |> filter(x < n_genes)

cum = ggplot(d, aes(pc, 100 * cum_ratio, colour = species)) +
  geom_vline(data = guides_df, aes(xintercept = x), linetype = 'dotted', colour = 'grey55') +
  geom_line(linewidth = 0.9) +
  facet_wrap(~tissue_lab, nrow = 1, scales = 'free_x') +
  scale_colour_manual(values = species_cols, name = 'Species') +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(breaks = \(l) { b = pretty(l); b[b >= 0 & b <= l[2] - 0.08 * diff(l)] }, expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = 'Number of principal components', y = 'Cumulative variance explained (%)') + theme_pc

scree = ggplot(filter(d, pc <= scree_pcs), aes(pc, 100 * var_ratio, colour = species)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  facet_wrap(~tissue_lab, nrow = 1) +
  scale_colour_manual(values = species_cols, name = 'Species', guide = 'none') +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = 'Principal component', y = 'Variance explained (%)') + theme_pc +
  theme(strip.text = element_blank())

p = (cum / scree) + plot_layout(heights = c(1.6, 1), guides = 'collect') +
  plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom', plot.tag = element_text(face = 'bold', size = 14))
ggsave(paste0(out, '.png'), p, width = 14, height = 7.5, dpi = 300, bg = 'white')
ggsave(paste0(out, '.pdf'), p, width = 14, height = 7.5, device = cairo_pdf, bg = 'white')

fmt = function(s) sprintf('%s %s: %.0f%% (10 PCs), %.0f%% (20), %s', tissue_labels[s$tissue], s$species, s$first10, s$first20,
                          if (s$n_genes <= 50) sprintf('100%% (all %d)', s$n_genes) else sprintf('%.0f%% (50)', s$first50))
writeLines(paste(
  'Explained variance of a separate PCA for each species within each tissue (dissected mouse data, evaluation cells).',
  'Each species\' log-normalised expression was scaled within the species (unit variance, clipped at 10) and all components kept,',
  'so each curve reaches 100% at the number of genes in the panel. A: cumulative variance explained; dotted lines at 10, 20 and 50 PCs.',
  sprintf('B: variance explained by each of the first %d PCs.', scree_pcs),
  paste0('Variance in the first 10 / 20 / 50 PCs: ', paste(sapply(seq_len(nrow(summ)), \(i) fmt(summ[i, ])), collapse = '; '), '.'),
  sep = '\n'), paste0(out, '_caption.txt'))
print(as.data.frame(summ |> mutate(across(c(PC1, first10, first20, first50), \(x) round(x, 1)))))
