# Selected metrics for the basic methods in each tissue (dissected):  Rscript analysis/basic_metrics.R
# Basic methods = every algorithm on its native input + LIGER UINMF; Unintegrated (joint log-normalised PCA) for reference.
# Raw metric values (all 0-1, higher = better): LISI-CoM (F1 of iLISI and cLISI, analysis/extra_metrics.py), ARI (scib, cluster vs
# label), SAS (analysis/sas_scores.py), iLISI (scib; identical to the extra-metrics iLISI), MAP (analysis/extra_metrics.py).
# Graph-only methods (BBKNN, SAMap) have no SAS or MAP: shown as dashes.
# Rows ordered as figures/literature/literature_heatmap (mean rank over tissues, scib-only overall improvement).
# Output: figures/literature/basic_metrics.{pdf,png,csv}
source('library.R')

# ---------------- config ----------------
out        = 'figures/literature'
extra_file = 'extra_metrics/extra_metrics_all.csv'
sas_file   = 'sas/sas_all.csv'
order_file = 'figures/literature/literature_comparison.csv'   # written by analysis/literature_comparison.R
metrics    = c(LISI_CoM = 'LISI-CoM', ARI = 'ARI', SAS = 'SAS', iLISI = 'iLISI', MAP = 'MAP')
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
tissue_short  = c(bone_marrow = 'BM', colon = 'Co', breast_cancer = 'BC', fetal_liver = 'FL')
highlight     = c('Harmony', 'Scanorama', 'CCA', 'RPCA')   # highlighted rows in the dot plot (Seurat = CCA, RPCA)
# -----------------------------------------
dir.create(out, showWarnings = FALSE, recursive = TRUE)
score_dir = function(t) if (t == 'breast_cancer') 'breast_cancer/data/dissected_mouse/' else sprintf('%s/data/05_evaluated_score/dissected_mouse/', t)
read_metrics = function(f) { d = suppressMessages(read_csv(f, show_col_types = FALSE)); tibble(metric = d[[1]], value = as.numeric(d[[2]]), method = names(d)[2]) }

ord = read_csv(order_file, show_col_types = FALSE)$method
rows = c(ord, 'Unintegrated')
scib = map_dfr(names(tissue_labels), \(t) map_dfr(list.files(score_dir(t), '^05\\..*\\.csv$', full.names = TRUE), read_metrics) |> mutate(tissue = t)) |>
  filter(metric %in% c('ARI_cluster/label', 'iLISI')) |> mutate(metric = recode(metric, `ARI_cluster/label` = 'ARI'))
extra = read_csv(extra_file, show_col_types = FALSE) |> filter(dataset %in% c('dissected', 'dissected_mouse')) |>
  select(tissue, method, LISI_CoM, MAP) |> pivot_longer(c(LISI_CoM, MAP), names_to = 'metric')
sas = read_csv(sas_file, show_col_types = FALSE) |> filter(dataset %in% c('dissected', 'dissected_mouse')) |>
  transmute(tissue, method, metric = 'SAS', value = sas)
d = bind_rows(scib, extra, sas) |> filter(method %in% rows) |> select(tissue, method, metric, value) |>
  complete(tissue = names(tissue_labels), method = rows, metric = names(metrics))
stopifnot(!any(duplicated(select(d, tissue, method, metric))))
# scVI-latent versions of the highlighted baselines (dot plot only)
dl = bind_rows(scib, extra, sas) |> filter(method %in% paste(highlight, 'scVI-latent')) |> select(tissue, method, metric, value)
write_csv(d |> arrange(metric, method, tissue), file.path(out, 'basic_metrics.csv'))
cat('Missing values (expected: SAS and MAP for BBKNN and SAMap):\n')
d |> filter(is.na(value)) |> dplyr::count(method, metric) |> as.data.frame() |> print(row.names = FALSE)

lab = function(m) ifelse(m %in% c('CCA', 'RPCA'), paste0(m, ' (Seurat)'), m)
yl = setNames(c(rev(seq_along(ord)) + 0.4, 0), rows)   # gap above the Unintegrated reference row
p = d |> mutate(metric = factor(metrics[metric], levels = metrics), x = factor(tissue_short[tissue], levels = tissue_short), y = yl[method]) |>
  ggplot(aes(x, y)) +
  geom_tile(aes(fill = value), colour = 'white', linewidth = 0.6, width = 0.96, height = 0.92) +
  geom_text(aes(label = ifelse(is.na(value), '–', sprintf('%.2f', value)), colour = !is.na(value) & value > 0.55), size = 2.7) +
  facet_wrap(~ metric, nrow = 1) +
  scale_fill_gradient(low = '#f2f0f7', high = '#3f007d', limits = c(0, 1), na.value = 'grey92', name = 'Value') +
  scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = 'grey15'), guide = 'none') +
  scale_x_discrete(position = 'top') +
  scale_y_continuous(breaks = yl, labels = lab(names(yl)), expand = expansion(add = 0.1)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = 'bold', size = 11), axis.text = element_text(colour = 'grey10'),
        panel.spacing.x = unit(0.8, 'lines'), legend.position = 'bottom', strip.placement = 'outside')
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('basic_metrics.', ext)), p, width = 12, height = 6, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')

# ---- dot-plot version: one dot per tissue, diamond = mean over tissues, ticks = Unintegrated in each tissue ----
tissue_cols   = c(bone_marrow = '#E69F00', colon = '#009E73', breast_cancer = '#CC79A7', fetal_liver = '#56B4E9')   # as analysis/benchmark_figures.R
tissue_shapes = c(bone_marrow = 16, colon = 17, breast_cancer = 15, fetal_liver = 18)
yd = setNames(rev(seq_along(ord)), ord)
dd = d |> filter(method %in% ord) |> mutate(metric = factor(metrics[metric], levels = metrics), y = yd[method],
                                            tissue = factor(tissue, levels = names(tissue_labels)))
ref = d |> filter(method == 'Unintegrated', !is.na(value)) |> mutate(metric = factor(metrics[metric], levels = metrics),
                                                                     tissue = factor(tissue, levels = names(tissue_labels)))
mn = dd |> group_by(metric, method, y) |> summarise(value = mean(value), .groups = 'drop') |> filter(!is.na(value))
na_rows = dd |> group_by(metric, method, y) |> filter(all(is.na(value))) |> distinct(metric, method, y)
best = mn |> group_by(metric) |> slice_max(value, n = 1) |> ungroup()   # highest mean over tissues among the basic methods
lat = dl |> group_by(method, metric) |> summarise(value = mean(value), .groups = 'drop') |>
  mutate(base = str_remove(method, ' scVI-latent'), y = yd[base] - 0.24, metric = factor(metrics[metric], levels = metrics))
cat('\nBest mean per metric (basic methods):\n'); print(as.data.frame(best |> mutate(value = round(value, 3))), row.names = FALSE)
cat('\nscVI-latent version minus native, mean over tissues:\n')
lat |> left_join(mn |> select(metric, base = method, native = value), by = c('metric', 'base')) |>
  transmute(base, metric, native = round(native, 3), scvi_latent = round(value, 3), diff = round(value - native, 3)) |>
  arrange(metric, base) |> as.data.frame() |> print(row.names = FALSE)
p2 = ggplot(dd, aes(value, y)) +
  geom_hline(yintercept = yd, colour = 'grey92', linewidth = 0.4) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = yd[highlight] - 0.42, ymax = yd[highlight] + 0.42, fill = '#f5e08a', alpha = 0.55) +
  geom_segment(data = ref, aes(x = value, xend = value, y = 0.4, yend = length(ord) + 0.6, colour = tissue),
               linewidth = 0.5, linetype = 'dashed', alpha = 0.7, inherit.aes = FALSE) +
  geom_point(aes(colour = tissue, shape = tissue), size = 2.2, alpha = 0.9, na.rm = TRUE) +
  geom_point(data = mn, aes(fill = 'Mean over tissues'), shape = 23, size = 2.8, colour = 'black') +
  geom_point(data = lat, aes(fill = 'Same algorithm on scVI-latent input\n(mean over tissues)'), shape = 23, size = 2.4, colour = '#7B5EA7') +
  geom_point(data = best, aes(y = y + 0.32, alpha = 'Highest mean in this metric'), shape = 8, size = 2.2, colour = 'black', stroke = 0.8) +
  geom_text(data = na_rows, aes(x = 0.5, label = 'not applicable'), colour = 'grey60', size = 2.7, fontface = 'italic') +
  facet_wrap(~ metric, nrow = 1) +
  scale_colour_manual(values = tissue_cols, labels = tissue_labels, name = NULL) +
  scale_shape_manual(values = tissue_shapes, labels = tissue_labels, name = NULL) +
  scale_fill_manual(values = c(`Mean over tissues` = 'black', `Same algorithm on scVI-latent input\n(mean over tissues)` = '#cbbde2'), name = NULL,
                    guide = guide_legend(order = 2, override.aes = list(shape = 23, size = 2.8, colour = c('black', '#7B5EA7')))) +
  scale_alpha_manual(values = c(`Highest mean in this metric` = 1), name = NULL, guide = guide_legend(order = 3)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '1')) +
  scale_y_continuous(breaks = yd, labels = lab(names(yd)), expand = expansion(add = 0.6)) +
  labs(x = 'Metric value (higher = better); dashed lines: Unintegrated in each tissue', y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(fill = 'grey95', colour = 'grey70'),
        strip.text = element_text(face = 'bold'), legend.position = 'bottom', legend.box = 'vertical', legend.margin = margin(0, 0, 0, 0),
        axis.text.y = element_text(colour = 'grey10', face = ifelse(names(yd) %in% highlight, 'bold', 'plain')))
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('basic_metrics_dots.', ext)), p2, width = 12.5, height = 6, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')
