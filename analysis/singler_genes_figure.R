# Genes SingleR uses to recognise each cell type, same-species reference (analysis/check_singler_genes.R: SingleR retrained as in
# the annotation, de.method = 't'; per type, the union over pairwise comparisons of the genes picked as up-regulated in it).
# Breast cancer mouse left blank. Reads saved CSVs only. Usage: Rscript45 analysis/singler_genes_figure.R
setwd(Sys.getenv('PROJECT_ROOT', '.')); source('library.R')
suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr); library(ggplot2) })
tissue_order = c('colon', 'bone_marrow', 'fetal_liver', 'breast_cancer')
tissue_name = \(t) str_to_sentence(str_replace(t, '_', ' '))
sg = list.files('figures/cross_tissue/singler_genes', '\\.csv$', full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |> bind_rows()
heads = sg |> distinct(tissue, species, n_shared_genes) |>
  mutate(panel = sprintf('%s, %s (%d panel genes)', tissue_name(tissue), species, n_shared_genes))
levels = c(as.vector(t(outer(tissue_order, c('human', 'mouse'), paste))))
heads = heads |> bind_rows(tibble(tissue = 'breast_cancer', species = 'mouse', panel = 'Breast cancer, mouse')) |>
  arrange(match(paste(tissue, species), levels))
d = sg |> inner_join(select(heads, tissue, species, panel), by = c('tissue', 'species')) |>
  mutate(panel = factor(panel, levels = heads$panel), y = paste(type, panel, sep = '___'))
d = mutate(d, y = factor(y, levels = y[order(panel, n_singler_genes, decreasing = c(TRUE, FALSE), method = 'radix')]))
write_csv(select(d, tissue, species, type, n_singler_genes, n_shared_genes, genes), 'figures/cross_tissue/singler_genes_per_type.csv')
p = ggplot(d, aes(n_singler_genes, y)) +
  geom_col(width = 0.7, fill = '#D2691E') +
  geom_text(aes(label = n_singler_genes), hjust = -0.3, size = 2.7, colour = 'grey25') +
  facet_wrap(~ panel, ncol = 2, scales = 'free_y', drop = FALSE) +
  scale_y_discrete(labels = \(x) sub('___.*$', '', x)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = 'Genes SingleR uses to recognise the cell type', y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = 'bold', hjust = 0))
ggsave('figures/cross_tissue/singler_genes_per_type.png', p, width = 11, height = 12, dpi = 300)
ggsave('figures/cross_tissue/singler_genes_per_type.pdf', p, width = 11, height = 12, device = cairo_pdf)
writeLines(paste('Number of genes SingleR uses to recognise each cell type of the same-species reference: the genes it selects as',
                 'up-regulated in the type in any pairwise comparison with another type (de.method = t, default number per comparison),',
                 'from the panel genes shared with the reference (in brackets). SingleR retrained exactly as in the annotation; the retrained',
                 'model reproduces the saved SingleR labels for 100% of 3,000 sampled cells per dataset.',
                 'Breast cancer mouse left blank: its reference has no Cancer Epithelial type.', sep = '\n'),
           'figures/cross_tissue/singler_genes_per_type_caption.txt')
cat('done figure\n')
