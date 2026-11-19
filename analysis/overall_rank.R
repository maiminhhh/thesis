# Overall marking of the integration methods on dissected tissue:  Rscript analysis/overall_rank.R
# Uses the existing scores (overall = 0.6 bio + 0.4 batch, unchanged), averaged over the four tissues
# (each tissue counts equally), in two versions:
#   rank:  rank within each tissue, then the mean rank;
#   score: the mean raw score.
# Point = full-data mean; small letters = the value in each tissue; rows sorted by the full-data overall mean.
# Bar = 95% interval over the 50 paired bootstrap subsamples. Rank: percentiles of the subsample mean ranks.
# Score: subsampling raises the scores by ~0.01-0.02, so the bar is the spread of the subsample means around
# their median, placed on the full-data mean (full + [q2.5 - median, q97.5 - median]).
# Graph-only methods (bootstrap not valid) and the unintegrated inputs are left out.
# Score definition (config): extra_bio / extra_batch = MAP / SAS (as analysis/benchmark_figures.R: bio = mean(NMI, ARI, ASW_label,
# cLISI, MAP), batch = mean(PCR_batch, ASW_label/batch, iLISI, graph_conn, SAS)), recomputed on every subsample
# (bootstrap_extra/extra_long.csv; full data: extra_metrics_all.csv, sas_all.csv); character(0) = scib only (original).
# Output: figures/benchmark_summary/15_overall_{rank,score}<suffix>.{png,csv}, 16_overall_heatmap<suffix>.{png,csv}
source('library.R')

# ---------------- config ----------------
in_file   = 'bootstrap/summary/1_long.csv'
out       = 'figures/benchmark_summary'
scores    = c(overall = 'Overall (0.6 bio + 0.4 batch)', bio = 'Bio conservation', batch = 'Batch correction')
exclude   = 'BBKNN|SAMap|scIB|Unintegrated'
highlight = c('Harmony', 'Scanorama', 'CCA', 'RPCA')
tissue_letters = c(bone_marrow = 'M', colon = 'C', breast_cancer = 'B', fetal_liver = 'F')
extra_bio     = c('MAP')        # character(0) for both: scib metrics only (the original outputs, no suffix)
extra_batch   = c('SAS')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI')
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
extra_boot    = 'bootstrap_extra/extra_long.csv'
extra_full    = c(MAP = 'extra_metrics/extra_metrics_all.csv', SAS = 'sas/sas_all.csv')
# -----------------------------------------
suffix = if (length(c(extra_bio, extra_batch))) paste0('_', paste(c(extra_bio, extra_batch), collapse = '_')) else ''

m = read_csv(in_file, show_col_types = FALSE, col_types = cols(`repeat` = 'c')) |>
  filter(metric %in% c(bio_metrics, batch_metrics), !str_detect(method, exclude)) |> select(dataset, `repeat`, method, metric, score)
if (length(c(extra_bio, extra_batch))) {
  ex = c(extra_bio, extra_batch)
  boot_x = read_csv(extra_boot, show_col_types = FALSE) |> filter(!str_detect(method, exclude)) |>
    transmute(dataset = tissue, `repeat` = as.character(rep), method, MAP, SAS) |>
    pivot_longer(all_of(ex), names_to = 'metric', values_to = 'score')
  full_x = bind_rows(read_csv(extra_full[['MAP']], show_col_types = FALSE) |> transmute(tissue, dataset, method, metric = 'MAP', score = MAP),
                     read_csv(extra_full[['SAS']], show_col_types = FALSE) |> transmute(tissue, dataset, method, metric = 'SAS', score = sas)) |>
    filter(dataset %in% c('dissected', 'dissected_mouse'), metric %in% ex, !str_detect(method, exclude)) |>
    transmute(dataset = tissue, `repeat` = 'point', method, metric, score)
  m = bind_rows(m, boot_x, full_x)
}
d = m |> filter(!is.na(score)) |>
  mutate(group = if_else(metric %in% c(bio_metrics, extra_bio), 'bio', 'batch')) |>
  group_by(dataset, `repeat`, method, group) |> summarise(score = mean(score), n = n(), .groups = 'drop')
n_need = c(bio = length(c(bio_metrics, extra_bio)), batch = length(c(batch_metrics, extra_batch)))
d = semi_join(d, d |> group_by(dataset, `repeat`, method) |>   # keep method x subsample with every metric present
                filter(all(n == n_need[group]), n() == 2) |> ungroup(), by = c('dataset', 'repeat', 'method')) |>
  select(-n) |> pivot_wider(names_from = group, values_from = score) |> mutate(overall = 0.6 * bio + 0.4 * batch) |>
  pivot_longer(c(overall, bio, batch), names_to = 'metric', values_to = 'score') |>
  group_by(dataset, `repeat`, metric) |> mutate(rank = rank(-score)) |> ungroup()
cat(sprintf('score definition: bio = %d metrics, batch = %d metrics; methods %d; subsamples per tissue: %s\n',
            n_need[['bio']], n_need[['batch']], n_distinct(d$method),
            paste(d |> filter(`repeat` != 'point', metric == 'overall') |> group_by(dataset) |>
                    summarise(k = n_distinct(`repeat`)) |> with(sprintf('%s %d', dataset, k)), collapse = ', ')))
# check: full-data scores = the benchmark scores of analysis/benchmark_figures.R (same definition)
bf = read_csv(sprintf('figures/benchmark_summary/benchmark_improvement%s.csv', suffix), show_col_types = FALSE) |>
  select(dataset = tissue, method, overall, bio, batch) |> pivot_longer(c(overall, bio, batch), names_to = 'metric', values_to = 'bf')
chk = inner_join(filter(d, `repeat` == 'point'), bf, by = c('dataset', 'method', 'metric'))
cat(sprintf('check vs analysis/benchmark_figures.R (full data): %d of %d matched, max |difference| = %.2g\n',
            nrow(chk), nrow(filter(d, `repeat` == 'point')), max(abs(chk$score - chk$bf))))

for (what in c('rank', 'score')) {
  dv = mutate(d, v = if (what == 'rank') rank else score)
  per_tissue = filter(dv, `repeat` == 'point') |> select(tissue = dataset, method, metric, v)
  full = per_tissue |> group_by(method, metric) |>
    summarise(mean = mean(v), worst = if (what == 'rank') max(v) else min(v), .groups = 'drop')
  boot = filter(dv, `repeat` != 'point') |> group_by(method, metric, `repeat`) |>
    summarise(m = mean(v), .groups = 'drop') |> group_by(method, metric) |>
    summarise(boot_median = median(m), q025 = quantile(m, 0.025), q975 = quantile(m, 0.975), n = n(), .groups = 'drop')
  res = full_join(full, boot, by = c('method', 'metric')) |>
    mutate(lo = if (what == 'rank') q025 else mean + q025 - boot_median,
           hi = if (what == 'rank') q975 else mean + q975 - boot_median)

  ord = res |> filter(metric == 'overall') |> arrange(if (what == 'rank') mean else -mean) |> pull(method)
  lab = function(m) ifelse(m %in% c('CCA', 'RPCA'), paste0(m, ' (Seurat)'), m)
  prep = function(x) x |> mutate(
    input = case_when(str_detect(method, 'UINMF') ~ 'UINMF', str_detect(method, 'scVI-latent') ~ 'scVI latent',
                      TRUE ~ 'PCA / own'),
    input = factor(input, levels = c('PCA / own', 'UINMF', 'scVI latent')),
    y = factor(lab(method), levels = rev(lab(ord))),
    panel = factor(scores[metric], levels = scores))
  res = prep(res); per_tissue = prep(per_tissue) |> mutate(letter = tissue_letters[tissue])

  write_csv(res |> select(metric, method, mean, worst, boot_median, boot_q025 = q025, boot_q975 = q975, lo, hi, n) |>
              arrange(metric, if (what == 'rank') mean else -mean),
            file.path(out, sprintf('15_overall_%s%s.csv', what, suffix)))

  input_cols = c('PCA / own' = '#2a78d6', 'UINMF' = '#eb6834', 'scVI latent' = '#1baf7a')
  n_m = length(ord)
  x_scale = if (what == 'rank') scale_x_continuous(limits = c(1, n_m), breaks = c(1, 5, 10, 15, n_m)) else
    scale_x_continuous(breaks = scales::breaks_pretty(4))
  bar_text = if (what == 'rank') 'Bar: 95% interval of that mean rank over ' else
    'Bar: 95% spread of the mean score over '
  p = ggplot(res, aes(y = y)) +
    geom_blank(aes(x = mean)) +   # sets the row order before the shading layer
    geom_hline(data = filter(res, method %in% highlight), aes(yintercept = y), colour = '#f5e08a', linewidth = 7, alpha = 0.6) +
    geom_text(data = per_tissue, aes(x = v, label = letter, colour = input), size = 2.6, alpha = 0.55, fontface = 'bold') +
    geom_errorbar(aes(xmin = lo, xmax = hi, colour = input), width = 0, orientation = 'y', linewidth = 1.1) +
    geom_point(aes(x = mean, fill = input), shape = 21, size = 3.2, colour = 'black', stroke = 0.5) +
    facet_wrap(~ panel, nrow = 1, scales = if (what == 'rank') 'fixed' else 'free_x') +
    x_scale +
    scale_fill_manual(values = input_cols, name = 'Input') +
    scale_colour_manual(values = input_cols, guide = 'none') +
    labs(x = if (what == 'rank') sprintf('Rank (1 = best of %d)', n_m) else 'Score (higher = better)', y = NULL,
         title = sprintf('Overall marking of the integration methods: mean %s over four tissues, dissected%s', what,
                         if (suffix != '') ' (+MAP +SAS)' else ''),
         caption = paste0('Point: mean of the four per-tissue ', what, 's on full data. ', bar_text,
                          max(res$n, na.rm = TRUE), ' paired subsamples of 20,000 cells (cell sampling only, not tissue choice)',
                          if (what == 'score') ',
placed on the full-data mean (subsampling raises the scores by ~0.01-0.02). ' else '.
',
                          'Letters: ', what, ' in each tissue (M bone marrow, C colon, B breast cancer, F fetal liver). ',
                          'Shaded: Harmony, Scanorama, Seurat (CCA, RPCA).\nRows sorted by overall mean ', what, '. ',
                          'Graph-only methods (BBKNN, SAMap) and unintegrated inputs left out.')) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.title.position = 'plot',
          strip.background = element_rect(fill = 'grey95', colour = 'grey70'), strip.text = element_text(face = 'bold'),
          legend.position = 'top', plot.caption = element_text(hjust = 0, colour = 'grey35'))
  ggsave(file.path(out, sprintf('15_overall_%s%s.png', what, suffix)), p, width = 15, height = 7.5, dpi = 300, bg = 'white')
}

# ---- main figure: heatmap table of the full-data scores ----
# Rows = methods sorted by mean overall score; columns = the four tissues and their mean, for overall, bio and batch.
# Number = the score; colour = its rank within the column (so tissues with different score ranges look alike).
# Output: figures/benchmark_summary/16_overall_heatmap.{png,csv}
tissue_names = c(bone_marrow = 'Bone\nmarrow', colon = 'Colon', breast_cancer = 'Breast\ncancer', fetal_liver = 'Fetal\nliver')
hm = filter(d, `repeat` == 'point') |> select(column = dataset, method, metric, score)
hm = bind_rows(hm, hm |> group_by(method, metric) |> summarise(column = 'mean', score = mean(score), .groups = 'drop')) |>
  group_by(metric, column) |> mutate(col_rank = rank(-score)) |> ungroup()
ord = hm |> filter(metric == 'overall', column == 'mean') |> arrange(-score) |> pull(method)
write_csv(hm |> select(metric, column, method, score, col_rank) |> arrange(metric, column, col_rank),
          file.path(out, sprintf('16_overall_heatmap%s.csv', suffix)))
lab = function(m) ifelse(m %in% c('CCA', 'RPCA'), paste0(m, ' (Seurat)'), m)
hm = hm |> mutate(
  x = factor(ifelse(column == 'mean', 'Mean', tissue_names[column]), levels = c(tissue_names, 'Mean')),
  y = factor(lab(method), levels = rev(lab(ord))),
  panel = factor(scores[metric], levels = scores),
  is_mean = column == 'mean')
n_m = length(ord)
hl_rows = filter(hm, method %in% highlight) |> distinct(y, panel)
lab_face = ifelse(levels(hm$y) %in% lab(highlight), 'bold', 'plain')
p = ggplot(hm, aes(x, y)) +
  geom_blank() +   # sets the row order before the shading layer
  geom_tile(data = hl_rows, aes(x = 3, y = y), width = Inf, height = 1, fill = '#f5e08a', alpha = 0.6, inherit.aes = FALSE) +
  geom_tile(aes(fill = col_rank), colour = 'white', linewidth = 0.6, width = 0.94, height = 0.9) +
  geom_tile(data = filter(hm, is_mean), fill = NA, colour = 'grey20', linewidth = 0.5, width = 0.94, height = 0.9) +
  geom_text(aes(label = sprintf('%.3f', score), colour = col_rank <= n_m * 0.45), size = 2.7,
            fontface = ifelse(hm$is_mean, 'bold', 'plain')) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_gradient(low = '#08306b', high = '#f2f7fd', limits = c(1, n_m), breaks = c(1, 5, 10, 15, n_m),
                      name = 'Rank in column\n(1 = best)') +
  scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = 'grey15'), guide = 'none') +
  scale_x_discrete(position = 'top') +
  labs(x = NULL, y = NULL,
       title = 'Integration methods on dissected tissue: scores per tissue and their mean',
       caption = paste0('Number: score on the full data (overall = 0.6 bio + 0.4 batch). Colour: rank within the column. ',
                        'Rows sorted by mean overall score.\nBold rows, shaded: Harmony, Scanorama, Seurat (CCA, RPCA). ',
                        'Graph-only methods (BBKNN, SAMap) and unintegrated inputs left out.\n',
                        'Uncertainty from cell sampling: see 15_overall_rank / 15_overall_score (bootstrap, 50 subsamples).')) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), plot.title.position = 'plot', strip.text = element_text(face = 'bold', size = 11),
        axis.text.y = element_text(face = lab_face, colour = 'grey10'), axis.text.x.top = element_text(colour = 'grey10'),
        legend.position = 'right', plot.caption = element_text(hjust = 0, colour = 'grey35'))
ggsave(file.path(out, sprintf('16_overall_heatmap%s.png', suffix)), p, width = 14, height = 7.5, dpi = 300, bg = 'white')
