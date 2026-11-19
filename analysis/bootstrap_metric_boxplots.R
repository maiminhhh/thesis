# Per-tissue bootstrap distributions of ARI, iLISI and the additional metrics:  Rscript analysis/bootstrap_metric_boxplots.R
# One figure per metric, one panel per tissue (dissected), methods sorted by bootstrap median.
# Box = middle 50% of the 50 paired subsamples (20k cells each); whiskers = 2.5th-97.5th percentiles (95% interval);
# diamond = the full-data score. Graph-only methods are left out (see below). Harmony, Scanorama and Seurat (CCA, RPCA) are shaded.
# Bootstrap values: overall, bio, batch (the saved scores, weights unchanged), ARI, iLISI from bootstrap/summary/1_long.csv; LISI-CoM and silhouette F1 per repeat from its
# iLISI/cLISI and ASW_label/batch / ASW_label (harmonic means, as analysis/extra_metrics.py); MAP, NC, SAS from
# bootstrap_extra/extra_long.csv (only metrics whose bootstrap exists are drawn).
# P(A > B) for every pair of the highlighted baselines: share of paired subsamples where A scores higher.
# Output: figures/benchmark_summary/14_bootstrap_<metric>.png, 14_bootstrap_intervals.csv, 14_bootstrap_baseline_pairs.csv
source('library.R')

out = 'figures/benchmark_summary'
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
dissected = c(bone_marrow = 'dissected', colon = 'dissected_mouse', breast_cancer = 'dissected_mouse', fetal_liver = 'dissected_mouse')
score_dir = c(bone_marrow = 'data/05_evaluated_score/dissected_mouse', colon = 'data/05_evaluated_score/dissected_mouse',
              breast_cancer = 'data/dissected_mouse', fetal_liver = 'data/05_evaluated_score/dissected_mouse')
graph_only = c('BBKNN', 'BBKNN UINMF', 'BBKNN scVI-latent', 'SAMap')
highlight = c('Harmony', 'Scanorama', 'CCA', 'RPCA')
metric_labels = c(overall = 'Overall score (0.6 bio + 0.4 batch)', bio = 'Bio conservation score',
                  batch = 'Batch correction score', ARI = 'ARI (bio)', MAP = 'MAP (bio)', NC = 'Neighbour conservation (bio)', iLISI = 'iLISI (batch)',
                  SAS = 'Seurat alignment score (batch)', LISI_CoM = 'LISI-CoM (combined)',
                  silhouette_F1 = 'Silhouette F1 (combined)')
f1 = function(a, b) 2 * a * b / (a + b)

# ---- bootstrap values ----
scib = read_csv('bootstrap/summary/1_long.csv', show_col_types = FALSE) |>
  filter(`repeat` != 'point', metric %in% c('ARI_cluster/label', 'iLISI', 'cLISI', 'ASW_label/batch', 'ASW_label', 'overall', 'bio', 'batch')) |>
  select(tissue = dataset, rep = `repeat`, method, metric, score) |>
  pivot_wider(names_from = metric, values_from = score) |>
  transmute(tissue, rep = as.integer(rep), method, overall, bio, batch, ARI = `ARI_cluster/label`, iLISI,
            LISI_CoM = f1(iLISI, cLISI), silhouette_F1 = f1(`ASW_label/batch`, ASW_label)) |>
  pivot_longer(c(overall, bio, batch, ARI, iLISI, LISI_CoM, silhouette_F1), names_to = 'metric', values_to = 'value')
boot = scib
if (file.exists('bootstrap_extra/extra_long.csv'))
  boot = bind_rows(boot, read_csv('bootstrap_extra/extra_long.csv', show_col_types = FALSE) |>
                     select(tissue, rep, method, MAP, NC, SAS) |>
                     pivot_longer(c(MAP, NC, SAS), names_to = 'metric', values_to = 'value'))
# graph-only methods are left out: restricting their stored graph to a 20k subsample (7-40% of the cells)
# breaks it (graph connectivity collapses, iLISI 0 or inflated), so their bootstrap values are not valid
boot = filter(boot, !is.na(value), !method %in% graph_only)

# ---- full-data values (as figure 13) ----
ari = imap_dfr(score_dir, function(d, t) map_dfr(Sys.glob(file.path(t, d, '05.*.csv')), function(f) {
  s = read.csv(f, row.names = 1, check.names = FALSE)
  tibble(tissue = t, method = names(s)[1], ARI = s['ARI_cluster/label', 1], iLISI = s['iLISI', 1])
}))
extra = read_csv('extra_metrics/extra_metrics_all.csv', show_col_types = FALSE) |>
  filter(dataset == dissected[tissue]) |> select(tissue, method, MAP, NC, LISI_CoM, silhouette_F1)
sas = read_csv('sas/sas_all.csv', show_col_types = FALSE) |> filter(dataset == dissected[tissue]) |> select(tissue, method, SAS = sas)
full = ari |> full_join(extra, by = c('tissue', 'method')) |> full_join(sas, by = c('tissue', 'method')) |>
  pivot_longer(-c(tissue, method), names_to = 'metric', values_to = 'full') |> filter(!is.na(full)) |>
  # overall, bio and batch on full data: the 'point' rows of 1_long (the saved scores)
  bind_rows(read_csv('bootstrap/summary/1_long.csv', show_col_types = FALSE) |>
              filter(`repeat` == 'point', metric %in% c('overall', 'bio', 'batch')) |>
              select(tissue = dataset, method, metric, full = score))

# ---- intervals ----
iv = boot |> group_by(tissue, method, metric) |>
  summarise(n = n(), q025 = quantile(value, 0.025), q25 = quantile(value, 0.25), median = median(value),
            q75 = quantile(value, 0.75), q975 = quantile(value, 0.975), .groups = 'drop') |>
  left_join(full, by = c('tissue', 'method', 'metric')) |>
  mutate(input = case_when(str_starts(method, 'Unintegrated') ~ 'Unintegrated',
                           str_detect(method, 'UINMF')        ~ 'UINMF',
                           str_detect(method, 'scVI-latent')  ~ 'scVI latent',
                           TRUE                               ~ 'PCA / own'),
         input = factor(input, levels = c('PCA / own', 'UINMF', 'scVI latent', 'Unintegrated')),
         label = case_when(method %in% graph_only ~ paste(method, '†'),
                           method %in% c('CCA', 'RPCA') ~ paste0(method, ' (Seurat)'), TRUE ~ method),
         panel = factor(tissue_labels[tissue], levels = tissue_labels),
         key = paste(label, tissue, sep = '___'))
write_csv(select(iv, metric, tissue, method, n, q025, q25, median, q75, q975, full), file.path(out, '14_bootstrap_intervals.csv'))

# ---- paired comparisons of the highlighted baselines ----
pairs = boot |> filter(method %in% highlight) |>
  pivot_wider(names_from = method, values_from = value) |>
  (\(w) map_dfr(combn(highlight, 2, simplify = FALSE), function(p) w |>
     filter(!is.na(.data[[p[1]]]), !is.na(.data[[p[2]]])) |>
     group_by(metric, tissue) |>
     summarise(A = p[1], B = p[2], n = n(), P_A_gt_B = mean(.data[[p[1]]] > .data[[p[2]]]),
               diff_median = median(.data[[p[1]]] - .data[[p[2]]]),
               diff_q025 = quantile(.data[[p[1]]] - .data[[p[2]]], 0.025),
               diff_q975 = quantile(.data[[p[1]]] - .data[[p[2]]], 0.975), .groups = 'drop')))() |>
  mutate(A = recode(A, CCA = 'CCA (Seurat)', RPCA = 'RPCA (Seurat)'), B = recode(B, CCA = 'CCA (Seurat)', RPCA = 'RPCA (Seurat)'),
         verdict = case_when(P_A_gt_B >= 0.95 ~ 'A > B', P_A_gt_B <= 0.05 ~ 'B > A', TRUE ~ 'unresolved'))
write_csv(pairs, file.path(out, '14_bootstrap_baseline_pairs.csv'))

# ---- one figure per metric ----
input_cols = c('PCA / own' = '#2a78d6', 'UINMF' = '#eb6834', 'scVI latent' = '#1baf7a', 'Unintegrated' = '#b5b5b0')
for (m in intersect(names(metric_labels), unique(iv$metric))) {
  x = iv |> filter(metric == m) |> mutate(key = fct_reorder(key, -median))
  hl = filter(x, method %in% highlight)
  n_rep = range(x$n)
  p = ggplot(x, aes(key)) +
    geom_blank(aes(y = median)) +   # sets the x order (by median) before the shading layer
    geom_vline(data = hl, aes(xintercept = key), colour = '#f5e08a', linewidth = 9, alpha = 0.6) +
    geom_boxplot(aes(ymin = q025, lower = q25, middle = median, upper = q75, ymax = q975, fill = input,
                     colour = method %in% highlight, linewidth = method %in% highlight),
                 stat = 'identity', width = 0.7) +
    geom_point(aes(y = full), shape = 23, size = 1.8, fill = 'white', colour = 'black', stroke = 0.6) +
    facet_wrap(~ panel, scales = 'free', ncol = 2) +
    scale_x_discrete(labels = function(k) str_remove(k, '___.*$')) +
    scale_fill_manual(values = input_cols, name = 'Input') +
    scale_colour_manual(values = c(`FALSE` = 'grey25', `TRUE` = 'black'), guide = 'none') +
    scale_linewidth_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.9), guide = 'none') +
    labs(x = NULL, y = metric_labels[[m]], title = sprintf('%s: bootstrap per tissue, dissected', metric_labels[[m]]),
         caption = paste0('Box: middle 50% of ', if (n_rep[1] == n_rep[2]) n_rep[1] else paste(n_rep, collapse = '-'),
                          ' paired subsamples of 20,000 cells; whiskers: 95% interval (2.5th-97.5th percentile); ',
                          'diamond: full-data score.\nShaded: Harmony, Scanorama and Seurat (CCA, RPCA). Sorted by ',
                          'bootstrap median, best on the left. Graph-only methods (BBKNN, SAMap) left out: their graph cannot be subsampled.\n',
                          'Overlapping boxes can still differ: the subsamples are paired. Subsampling shifts some metrics (e.g. iLISI) away from the full-data score.')) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), plot.title.position = 'plot',
          strip.background = element_rect(fill = 'grey95', colour = 'grey70'), strip.text = element_text(face = 'bold'),
          legend.position = 'top', plot.caption = element_text(hjust = 0, colour = 'grey35'),
          axis.text.x = element_text(size = 8.5, angle = 55, hjust = 1))
  ggsave(file.path(out, sprintf('14_bootstrap_%s.png', m)), p, width = 16, height = 11, dpi = 300, bg = 'white')
}
