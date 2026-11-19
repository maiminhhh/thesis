# Bar plot of cells per SingleR (cross-species reference) label in one dataset; counts from analysis/singler_cross_counts.py.
#   Rscript --vanilla analysis/singler_cross_barplot.R fetal_liver human
# Output: figures/<tissue>/annotation/singler_cross_counts_<species>.{png,pdf}
suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2)})
args = commandArgs(TRUE); tissue = args[1]; species = args[2]
species_cols = c(human = '#4477AA', mouse = '#EE6677')                     # as analysis/visualisations.R
stem = sprintf('figures/%s/annotation/singler_cross_counts_%s', tissue, species)
d = read_csv(paste0(stem, '.csv'), show_col_types = FALSE) |>
  mutate(kind = ifelse(label == 'Unassigned', 'Unassigned (pruned call)', 'Cell type'),
         label = factor(label, levels = rev(label[order(-n_cells, label)])),
         txt = ifelse(n_cells == 0, 'never assigned', sprintf('%s (%s%%)', scales::comma(n_cells), ifelse(pct > 0 & pct < 0.05, '<0.1', sprintf('%.1f', pct)))))
p = ggplot(d, aes(pct, label, fill = kind)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = txt), hjust = -0.08, size = 3, colour = 'grey20') +
  scale_fill_manual(values = c('Cell type' = species_cols[[species]], 'Unassigned (pruned call)' = 'grey65'), name = NULL) +
  scale_x_continuous(limits = c(0, max(d$pct) * 1.28), breaks = seq(0, 100, 10), labels = \(x) paste0(x, '%'), expand = c(0, 0)) +
  labs(x = sprintf('%% of cells (n = %s)', scales::comma(sum(d$n_cells))), y = NULL,
       title = sprintf('SingleR, cross-species reference: %s %s', tolower(sub('_', ' ', tissue)), species)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), legend.position = 'top',
        legend.justification = 'left', plot.title = element_text(face = 'bold', size = 12))
ggsave(paste0(stem, '.png'), p, width = 7.5, height = 5, dpi = 300, bg = 'white')
ggsave(paste0(stem, '.pdf'), p, width = 7.5, height = 5, bg = 'white', device = cairo_pdf)
cat('wrote', paste0(stem, '.{png,pdf}\n'))
