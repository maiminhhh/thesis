# Vertical boxplots of the extra metrics, ARI and iLISI, one panel per metric:  Rscript new_metrics_bars.R
# Dissected tissue, every method; boxplot over the four tissues (sorted by median), points = each tissue.
# ARI from the saved scib scores; iLISI (the saved scib value), MAP, NC, LISI-CoM, silhouette F1 from extra_metrics/extra_metrics_all.csv;
# SAS from sas/sas_all.csv. Graph-only methods (†) have no embedding: no MAP, NC, SAS or silhouette F1,
# and their LISI-CoM uses graph LISI. Higher is better for every metric.
# Output: figures/benchmark_summary/13_new_metrics_boxplots.{png,csv}
source('library.R')

out = 'figures/benchmark_summary'
dissected = c(bone_marrow = 'dissected', colon = 'dissected_mouse', breast_cancer = 'dissected_mouse', fetal_liver = 'dissected_mouse')
score_dir = c(bone_marrow = 'data/05_evaluated_score/dissected_mouse', colon = 'data/05_evaluated_score/dissected_mouse',
              breast_cancer = 'data/dissected_mouse', fetal_liver = 'data/05_evaluated_score/dissected_mouse')
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
graph_only = c('BBKNN', 'BBKNN UINMF', 'BBKNN scVI-latent', 'SAMap')
highlight = c('Harmony', 'Scanorama', 'CCA', 'RPCA')   # shaded; CCA and RPCA are Seurat's integrations

ari = imap_dfr(score_dir, function(d, t) map_dfr(Sys.glob(file.path(t, d, '05.*.csv')), function(f) {
  s = read.csv(f, row.names = 1, check.names = FALSE)
  tibble(tissue = t, method = names(s)[1], ARI = s['ARI_cluster/label', 1])
}))
extra = read_csv('extra_metrics/extra_metrics_all.csv', show_col_types = FALSE) |>
  filter(dataset == dissected[tissue]) |> select(tissue, method, iLISI, MAP, NC, LISI_CoM, silhouette_F1)
sas = read_csv('sas/sas_all.csv', show_col_types = FALSE) |>
  filter(dataset == dissected[tissue]) |> select(tissue, method, SAS = sas)

metric_labels = c(ARI = 'ARI (bio)', MAP = 'MAP (bio)', NC = 'Neighbour conservation (bio)',
                  iLISI = 'iLISI (batch)', SAS = 'Seurat alignment score (batch)', LISI_CoM = 'LISI-CoM (combined)',
                  silhouette_F1 = 'Silhouette F1 (combined)')
d = ari |> full_join(extra, by = c('tissue', 'method')) |> full_join(sas, by = c('tissue', 'method')) |>
  pivot_longer(names(metric_labels), names_to = 'metric', values_to = 'value') |>
  filter(!is.na(value)) |>
  mutate(input = case_when(str_starts(method, 'Unintegrated') ~ 'Unintegrated',
                           str_detect(method, 'UINMF')        ~ 'UINMF',
                           str_detect(method, 'scVI-latent')  ~ 'scVI latent',
                           TRUE                               ~ 'PCA / own'),
         input = factor(input, levels = c('PCA / own', 'UINMF', 'scVI latent', 'Unintegrated')),
         label = case_when(method %in% graph_only ~ paste(method, '†'),
                           method %in% c('CCA', 'RPCA') ~ paste0(method, ' (Seurat)'),
                           TRUE ~ method),
         tissue = factor(tissue_labels[tissue], levels = tissue_labels),
         metric = factor(metric_labels[metric], levels = metric_labels))
stopifnot(all(dplyr::count(d, metric, method)$n == 4))   # every method scored in all four tissues

m = d |> group_by(metric, method, label, input) |>
  summarise(median = median(value), mean = mean(value), .groups = 'drop') |>
  mutate(key = paste(label, metric, sep = '___'), key = fct_reorder(key, -median))
d = d |> mutate(key = factor(paste(label, metric, sep = '___'), levels = levels(m$key)))
write_csv(select(m, metric, method, input, median, mean), file.path(out, '13_new_metrics_boxplots.csv'))

input_cols = c('PCA / own' = '#2a78d6', 'UINMF' = '#eb6834', 'scVI latent' = '#1baf7a', 'Unintegrated' = '#b5b5b0')
hl = filter(m, method %in% highlight)
p = ggplot(d, aes(key, value)) +
  geom_blank() +   # sets the x order (by median) before the shading layer
  geom_vline(data = hl, aes(xintercept = key), colour = '#f5e08a', linewidth = 9, alpha = 0.6) +
  geom_boxplot(aes(fill = input, colour = method %in% highlight, linewidth = method %in% highlight),
               width = 0.7, outlier.shape = NA) +
  scale_colour_manual(values = c(`FALSE` = 'grey25', `TRUE` = 'black'), guide = 'none') +
  scale_linewidth_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.9), guide = 'none') +
  geom_point(aes(shape = tissue), size = 1.4, colour = 'grey10', alpha = 0.8) +
  facet_wrap(~ metric, scales = 'free', ncol = 2) +
  scale_x_discrete(labels = function(x) str_remove(x, '___.*$')) +
  scale_fill_manual(values = input_cols, name = 'Input') +
  scale_shape_manual(values = c(16, 17, 15, 18), name = 'Tissue') +
  labs(x = NULL, y = 'Score (higher is better)',
       title = 'Methods on ARI, iLISI and the additional metrics, dissected tissue',
       caption = paste('Box: the four tissues (sorted by median, best on the left); points: each tissue. Shaded: Harmony, Scanorama and Seurat (CCA, RPCA).\n † graph-only method: no embedding, so no',
                       'MAP, NC, SAS or silhouette F1; its iLISI and LISI-CoM use graph LISI.\nLISI-CoM and silhouette F1 are',
                       'harmonic means of the saved scib iLISI/cLISI and batch/cell-type ASW. NC: overlap with each',
                       "species' own PCA neighbours.")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), plot.title.position = 'plot',
        strip.background = element_rect(fill = 'grey95', colour = 'grey70'), strip.text = element_text(face = 'bold'),
        legend.position = 'top', plot.caption = element_text(hjust = 0, colour = 'grey35'),
        axis.text.x = element_text(size = 8.5, angle = 55, hjust = 1))
ggsave(file.path(out, '13_new_metrics_boxplots.png'), p, width = 18, height = 22, dpi = 300, bg = 'white')
