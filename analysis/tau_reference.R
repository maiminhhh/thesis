# Reference for the label-sensitivity taus: how much do method rankings change between datasets of the same
# tissue under the same (EBCC) labels? Kendall's tau between the EBCC-label rankings of every pair of datasets,
# for bio-conservation and overall (0.6 bio + 0.4 batch), next to the label-change taus
# (figures/cross_tissue/label_sensitivity.csv). Reads saved 05 scores only.
# Usage: Rscript analysis/tau_reference.R  -> figures/cross_tissue/tau_reference<suffix>.csv
# matched = TRUE: the same scores and methods as analysis/label_robustness.R (bio = NMI, ARI, ASW_label, cLISI; batch = PCR_batch,
# ASW_label/batch, iLISI, graph_conn; the 23 methods, unintegrated baselines excluded), also for batch correction, and
# compared with label_robustness.csv (colon excluded, as there); outputs tau_reference_matched*.csv.
# matched = FALSE: the original (bio with the two isolated-label metrics, unintegrated baselines included,
# compared with label_sensitivity.csv).
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(readr); library(stringr) })
setwd(Sys.getenv('PROJECT_ROOT', '.'))

# ---------------- config ----------------
matched = TRUE
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = if (matched) c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI') else
  c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'isolated_label_F1', 'isolated_label_silhouette', 'cLISI')
# matched only: extra = TRUE uses the primary +MAP +SAS score (bio + MAP, batch + SAS; graph-only methods without them, as
# analysis/benchmark_figures.R) on the datasets with MAP / SAS (no full tissue), compared with label_robustness_MAP_SAS.csv
extra = TRUE
suffix = paste0(if (matched) '_matched' else '', if (matched && extra) '_MAP_SAS' else '')
# -----------------------------------------
extra = matched && extra
datasets = c(if (!extra) 'full', 'dissected', 'genes2', 'genes4', 'genes6', 'genes8')
lab_file = sprintf('figures/benchmark_summary/label_robustness%s.csv', if (extra) '_MAP_SAS' else '')
extra_long = if (!extra) tibble() else bind_rows(
  read_csv('extra_metrics/extra_metrics_all.csv', show_col_types = FALSE) |> transmute(tissue, dataset, method, metric = 'MAP', value = MAP),
  read_csv('sas/sas_all.csv', show_col_types = FALSE) |> transmute(tissue, dataset, method, metric = 'SAS', value = sas)) |>
  mutate(dataset = str_remove(dataset, '_mouse$')) |> filter(!is.na(value))
score_dir = function(tissue, d) {
  dn = if (d %in% c('full', 'dissected')) paste0(d, '_mouse') else d
  if (tissue == 'breast_cancer' || startsWith(d, 'genes')) sprintf('%s/data/%s/', tissue, dn)
  else sprintf('%s/data/05_evaluated_score/%s/', tissue, dn)
}
set_scores = function(dir, extra_m = NULL) {
  files = list.files(dir, pattern = '^05\\..*\\.csv$', full.names = TRUE)
  map_dfr(files, function(f) {
    df = read_csv(f, show_col_types = FALSE)
    tibble(metric = df[[1]], value = as.numeric(df[[2]]), method = names(df)[2])
  }) |> bind_rows(extra_m) |>
    mutate(group = case_when(metric %in% c(batch_metrics, if (extra) 'SAS') ~ 'batch',
                             metric %in% c(bio_metrics, if (extra) 'MAP') ~ 'bio')) |>
    filter(!is.na(group), !is.na(value)) |>
    group_by(method, group) |> summarise(score = mean(value), .groups = 'drop') |>
    pivot_wider(names_from = group, values_from = score) |>
    mutate(overall = 0.6 * bio + 0.4 * batch)
}
tissues = c('colon', 'bone_marrow', 'fetal_liver', 'breast_cancer')
scores = map_dfr(tissues, \(t) map_dfr(datasets, \(d) mutate(set_scores(score_dir(t, d),
  if (extra) filter(extra_long, tissue == t, dataset == d) |> select(metric, value, method)), tissue = t, dataset = d)))
if (extra) {
  scores = filter(scores, !str_starts(method, 'Unintegrated'))
  # check: the dissected scores equal analysis/benchmark_figures.R's (+MAP +SAS)
  main = read_csv('figures/benchmark_summary/benchmark_improvement_MAP_SAS.csv', show_col_types = FALSE)
  chk = inner_join(filter(scores, dataset == 'dissected'), main, by = c('tissue', 'method'))
  stopifnot(nrow(chk) == 4 * 23, max(abs(chk$overall.x - chk$overall.y)) < 1e-12)
} else if (matched) {
  scores = filter(scores, !str_starts(method, 'Unintegrated'))
  # check: the EBCC scores equal those of analysis/label_robustness.R / benchmark_summary (scores_all.csv)
  main = read_csv('figures/benchmark_summary/scores_all.csv', show_col_types = FALSE)
  chk = inner_join(scores, main, by = c('tissue', 'dataset' = 'case', 'method'))
  stopifnot(nrow(chk) == 24 * 23, max(abs(chk$overall.x - chk$overall.y)) < 1e-12)
}
# Weighted Kendall tau (Vigna 2015), as analysis/mcb_analysis.R / analysis/benchmark_figures.R: hyperbolic weight 1/(r + 1) at importance
# rank r (0 = best), pair weight = sum of the two, averaged over both rankings as the importance order (= scipy weightedtau)
wtau = function(x, y) {
  n = length(x)
  one = function(a, b) {
    rho = integer(n); rho[order(-a, -b)] = 0:(n - 1)
    w = 1 / (rho + 1); W = outer(w, w, '+'); up = upper.tri(W)
    sa = sign(outer(a, a, '-')); sb = sign(outer(b, b, '-'))
    sum((W * sa * sb)[up]) / sqrt(sum((W * (sa != 0))[up]) * sum((W * (sb != 0))[up]))
  }
  (one(x, y) + one(y, x)) / 2
}
pairs = combn(datasets, 2, simplify = FALSE)
ref = map_dfr(tissues, \(t) map_dfr(pairs, function(p) {
  m = inner_join(filter(scores, tissue == t, dataset == p[1]), filter(scores, tissue == t, dataset == p[2]),
                 by = 'method', suffix = c('_a', '_b'))
  tibble(tissue = t, a = p[1], b = p[2], n_methods = nrow(m),
         adjacent = match(p[2], datasets) - match(p[1], datasets) == 1,
         tau_bio = cor(m$bio_a, m$bio_b, method = 'kendall'),
         tau_overall = cor(m$overall_a, m$overall_b, method = 'kendall'),
         tau_batch = cor(m$batch_a, m$batch_b, method = 'kendall'),
         wtau_bio = wtau(m$bio_a, m$bio_b), wtau_overall = wtau(m$overall_a, m$overall_b), wtau_batch = wtau(m$batch_a, m$batch_b))
}))
if (!matched) ref = select(ref, -tau_batch, -wtau_batch)
write_csv(ref, sprintf('figures/cross_tissue/tau_reference%s.csv', suffix))

lab = if (matched) read_csv(lab_file, show_col_types = FALSE) |>
  group_by(tissue, score) |> summarise(label_min = min(tau), label_median = median(tau), label_max = max(tau), .groups = 'drop') else
  read_csv('figures/cross_tissue/label_sensitivity.csv', show_col_types = FALSE) |>
  mutate(score = if_else(str_starts(score, 'Bio'), 'bio', 'overall')) |>
  group_by(tissue, score) |> summarise(label_min = min(tau), label_median = median(tau), label_max = max(tau), .groups = 'drop')
summ = ref |> pivot_longer(starts_with('tau_'), names_to = 'score', values_to = 'tau', names_prefix = 'tau_') |>
  group_by(tissue, score) |>
  summarise(n_methods = paste(range(n_methods), collapse = '-'),
            between_all_min = min(tau), between_all_median = median(tau),
            between_adjacent_median = median(tau[adjacent]), between_adjacent_max = max(tau[adjacent]), .groups = 'drop') |>
  left_join(lab, by = c('tissue', 'score')) |>
  mutate(tissue = factor(tissue, levels = tissues)) |> arrange(score, tissue)
write_csv(summ, sprintf('figures/cross_tissue/tau_reference%s_summary.csv', suffix))
print(as.data.frame(summ), digits = 2)
if (matched) {   # pooled over the three tissues with an alternative labelling (colon excluded), as in analysis/label_robustness.R
  # statistic: Kendall's tau ('tau') and weighted Kendall's tau ('weighted_tau'), for both comparisons
  pooled = bind_rows(
    ref |> filter(tissue != 'colon') |> pivot_longer(matches('^w?tau_'), names_to = c('statistic', 'score'), names_sep = '_',
                                                     values_to = 'tau') |>
      mutate(statistic = recode(statistic, wtau = 'weighted_tau'), comparison = 'between datasets (EBCC)'),
    read_csv(lab_file, show_col_types = FALSE) |> select(score, tau, weighted_tau) |>
      pivot_longer(c(tau, weighted_tau), names_to = 'statistic', values_to = 'tau') |>
      mutate(comparison = 'EBCC vs best annotator')) |>
    group_by(statistic, comparison, score) |>
    summarise(n = n(), min = min(tau), q25 = quantile(tau, 0.25), median = median(tau), q75 = quantile(tau, 0.75), max = max(tau), .groups = 'drop')
  write_csv(pooled, sprintf('figures/cross_tissue/tau_reference%s_pooled.csv', suffix))
  cat(sprintf('\nPooled over bone marrow, breast cancer and fetal liver (23 methods, %s score):\n', if (extra) '+MAP +SAS' else 'scib-only'))
  print(as.data.frame(pooled), digits = 2)
}
