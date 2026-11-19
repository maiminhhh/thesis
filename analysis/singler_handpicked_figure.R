# Hand-picked markers that SingleR also finds useful: per cell type, the hand-picked markers measured by the panel, split by
# whether SingleR uses them to identify that type, uses them only for other types, or does not use them
# (analysis/check_singler_handpicked.R -> figures/cross_tissue/singler_handpicked.csv). Same-species reference; breast cancer mouse
# left blank. Reads saved CSVs only. Usage: Rscript45 analysis/singler_handpicked_figure.R
setwd(Sys.getenv('PROJECT_ROOT', '.')); source('library.R')
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(stringr); library(ggplot2) })
tissue_order = c('colon', 'bone_marrow', 'fetal_liver', 'breast_cancer')
tissue_name = \(t) str_to_sentence(str_replace(t, '_', ' '))
d = read_csv('figures/cross_tissue/singler_handpicked.csv', show_col_types = FALSE)
heads = d |> group_by(tissue, species) |>
  summarise(used = sum(n_in_type_set), n = sum(n_handpicked_in_panel), .groups = 'drop') |>
  mutate(panel = sprintf('%s, %s: %d of %d used for their type (%.0f%%)', tissue_name(tissue), species, used, n, 100 * used / n)) |>
  bind_rows(tibble(tissue = 'breast_cancer', species = 'mouse', panel = 'Breast cancer, mouse')) |>
  arrange(match(tissue, tissue_order), species)
status_lab = c(this = 'Used by SingleR for this type', other = 'Used by SingleR only for other types', none = 'Not used by SingleR')
d = d |> inner_join(select(heads, tissue, species, panel), by = c('tissue', 'species')) |>
  mutate(panel = factor(panel, levels = heads$panel), this = n_in_type_set, other = n_in_shared_set - n_in_type_set,
         none = n_handpicked_in_panel - n_in_shared_set, y = paste(type, panel, sep = '___'),
         label = sprintf('%d/%d', n_in_type_set, n_handpicked_in_panel))
d = mutate(d, y = factor(y, levels = y[order(panel, n_handpicked_in_panel, n_in_type_set, decreasing = c(TRUE, FALSE, FALSE), method = 'radix')]))
long = d |> pivot_longer(c(this, other, none), names_to = 'status', values_to = 'n') |>
  mutate(status = factor(status_lab[status], levels = rev(unname(status_lab))))
p = ggplot(long, aes(n, y, fill = status)) +
  geom_col(width = 0.72, colour = 'white', linewidth = 0.25) +
  geom_text(data = d, aes(x = n_handpicked_in_panel, y = y, label = label), inherit.aes = FALSE, hjust = -0.2, size = 2.6, colour = 'grey25') +
  facet_wrap(~ panel, ncol = 2, scales = 'free', drop = FALSE) +
  scale_y_discrete(labels = \(x) sub('___.*$', '', x)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)), breaks = scales::breaks_pretty(4)) +
  scale_fill_manual(values = setNames(c('#D2691E', '#F2C18D', 'grey85'), unname(status_lab)), breaks = unname(status_lab), name = NULL) +
  labs(x = 'Hand-picked markers measured by the panel', y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = 'top', panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = 'bold', hjust = 0))
ggsave('figures/cross_tissue/singler_handpicked.png', p, width = 12, height = 13, dpi = 300)
ggsave('figures/cross_tissue/singler_handpicked.pdf', p, width = 12, height = 13, device = cairo_pdf)
writeLines(paste('Hand-picked markers of each cell type (the lists used for marker AUROC) that the gene panel measures, split by whether SingleR',
                 '(same-species reference, retrained exactly as in the annotation) uses them to identify that type: SingleR picks the 10 genes',
                 'most up-regulated in the type against each other type; a marker it picks only for other types still enters the correlation',
                 'over all types\' genes. Labels: used for this type / hand-picked markers in the panel. Only types in the SingleR reference.',
                 'Breast cancer mouse left blank: its reference has no Cancer Epithelial type.', sep = '\n'),
           'figures/cross_tissue/singler_handpicked_caption.txt')
cat('done figure\n')
