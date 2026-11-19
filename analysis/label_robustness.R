# Rank robustness to the labelling:  Rscript analysis/label_robustness.R
# Two labellings per dataset: EBCC (main benchmark) and the best single annotator
# (per species, highest mean marker AUROC; label_sensitivity_choice.json), scored by
# the same 05 metrics in <scores>/label_best_annotator/. Colon has no alternative:
# EBCC is the best annotator in both species, so its rankings are identical by construction.
# Scores use the benchmark definitions (analysis/benchmark_summary.R): 4 bio, 4 batch metrics,
# overall = 0.6 bio + 0.4 batch, 23 methods (Unintegrated baselines excluded).
# Output: figures/benchmark_summary/12_label_robustness*.{png,pdf} (overall ranking only; caption in *_caption.txt) and label_robustness*.csv (all three scores)
suppressMessages({library(tidyverse); library(ggrepel); library(patchwork)})

# ---------------- config ----------------
# Score definition. character(0) for both: the scib-only score (4 bio, 4 batch; the original outputs, no suffix, all
# six datasets). extra_bio = 'MAP', extra_batch = 'SAS': the primary +MAP +SAS score (as analysis/benchmark_figures.R: bio and
# batch = equal-weight means of the metrics that exist for the method; graph-only methods have no MAP / SAS), on the
# datasets with MAP / SAS (dissected, genes2-8; no full tissue). EBCC MAP: analysis/extra_metrics.py; best-annotator MAP:
# analysis/label_map_alt.py (same neighbours, other labels); SAS is label-free, so the same value under both labellings.
extra_bio   = c('MAP')
extra_batch = c('SAS')
map_files   = c(ebcc = 'extra_metrics/extra_metrics_all.csv', best_annotator = 'label_map_alt/map_alt_all.csv')
sas_file    = 'sas/sas_all.csv'
fig_stat    = 'weighted_tau'   # statistic in the figure and caption: 'weighted_tau' (as section 2.2) or 'tau'
ref_file    = 'figures/cross_tissue/tau_reference_matched_MAP_SAS_pooled.csv'   # analysis/tau_reference.R: tau between datasets of a
                                                                                 # tissue (EBCC), drawn behind panel A; skipped if absent
# -----------------------------------------
extra  = length(c(extra_bio, extra_batch)) > 0
suffix = if (extra) paste0('_', paste(c(extra_bio, extra_batch), collapse = '_')) else ''

out = 'figures/benchmark_summary'
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon',
                  breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
case_levels = c(if (!extra) 'full', 'dissected', 'genes2', 'genes4', 'genes6', 'genes8')
score_dir = function(tissue, case) {
  if (tissue == 'breast_cancer')
    return(sprintf('breast_cancer/data/%s/', if (case %in% c('full', 'dissected')) paste0(case, '_mouse') else case))
  if (case %in% c('full', 'dissected')) sprintf('%s/data/05_evaluated_score/%s_mouse/', tissue, case)
  else sprintf('%s/data/%s/', tissue, case)
}
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI')
read_metrics = function(f) {
  df = suppressMessages(read_csv(f, show_col_types = FALSE))
  tibble(metric = df[[1]], value = as.numeric(df[[2]]), method = names(df)[2])
}
# MAP / SAS per labelling, tissue, case, method (long; empty without extra metrics)
extra_long = if (!extra) tibble() else {
  sas = read_csv(sas_file, show_col_types = FALSE) |> transmute(tissue, case = str_remove(dataset, '_mouse$'), method, metric = 'SAS', value = sas)
  bind_rows(lapply(names(map_files), function(l) {
    mp = read_csv(map_files[[l]], show_col_types = FALSE)
    mp = if (l == 'ebcc') transmute(mp, tissue, case = str_remove(dataset, '_mouse$'), method, metric = 'MAP', value = MAP) else
      transmute(mp, tissue, case = str_remove(dataset, '_mouse$'), method, metric = 'MAP', value = MAP_best_annotator)
    bind_rows(mp, sas) |> filter(metric %in% c(extra_bio, extra_batch), !is.na(value)) |> mutate(labelling = l)
  }))
}
set_scores = function(dir, extra_m = NULL) {
  files = list.files(dir, pattern = '^05\\..*\\.csv$', full.names = TRUE)
  if (length(files) == 0) return(NULL)
  map_dfr(files, read_metrics) |> bind_rows(extra_m) |>
    mutate(group = case_when(metric %in% c(batch_metrics, extra_batch) ~ 'batch',
                             metric %in% c(bio_metrics, extra_bio)     ~ 'bio')) |>
    filter(!is.na(group), !is.na(value)) |>
    group_by(method, group) |> summarise(score = mean(value), n_metrics = n(), .groups = 'drop') |>
    pivot_wider(names_from = group, values_from = c(score, n_metrics)) |>
    rename(bio = score_bio, batch = score_batch) |>
    mutate(overall = 0.6 * bio + 0.4 * batch)
}
ex_for = function(l, t, cs) if (!extra) NULL else filter(extra_long, labelling == l, tissue == t, case == cs) |> select(metric, value, method)

scores = expand_grid(tissue = names(tissue_labels), case = case_levels) |>
  pmap_dfr(function(tissue, case) {
    d = score_dir(tissue, case)
    alt_dir = file.path(d, 'label_best_annotator')
    bind_rows(ebcc = set_scores(d, ex_for('ebcc', tissue, case)),
              best_annotator = if (dir.exists(alt_dir)) set_scores(alt_dir, ex_for('best_annotator', tissue, case)),
              .id = 'labelling') |> mutate(tissue = tissue, case = case)
  }) |>
  filter(!str_starts(method, 'Unintegrated'))

# the EBCC scores must reproduce the main benchmark exactly
if (!extra) {
  main = read_csv(file.path(out, 'scores_all.csv'), show_col_types = FALSE)
  chk = scores |> filter(labelling == 'ebcc') |> inner_join(main, by = c('tissue', 'case', 'method'))
  stopifnot(nrow(chk) == 24 * 23, max(abs(chk$overall.x - chk$overall.y)) < 1e-12)
} else {
  # dissected: benchmark_improvement<suffix>.csv (analysis/benchmark_figures.R); every embedding method needs MAP and SAS under
  # both labellings (the n_metrics of the two labellings must agree)
  main = read_csv(sprintf('%s/benchmark_improvement%s.csv', out, suffix), show_col_types = FALSE)
  chk = scores |> filter(labelling == 'ebcc', case == 'dissected') |> inner_join(main, by = c('tissue', 'method'))
  stopifnot(nrow(chk) == 4 * 23, max(abs(chk$overall.x - chk$overall.y)) < 1e-12,
            all(chk$n_metrics_bio.x == chk$n_metrics_bio.y), all(chk$n_metrics_batch.x == chk$n_metrics_batch.y))
  nm = scores |> select(labelling, tissue, case, method, n_metrics_bio, n_metrics_batch) |>
    pivot_wider(names_from = labelling, values_from = c(n_metrics_bio, n_metrics_batch)) |> filter(!is.na(n_metrics_bio_best_annotator))
  stopifnot(all(nm$n_metrics_bio_ebcc == nm$n_metrics_bio_best_annotator), all(nm$n_metrics_batch_ebcc == nm$n_metrics_batch_best_annotator),
            all(nm$n_metrics_bio_ebcc[!str_detect(nm$method, 'BBKNN|SAMap')] == length(bio_metrics) + length(extra_bio)))
  cat(sprintf('check: dissected EBCC scores = benchmark_improvement%s.csv (%d rows, max |diff| %.2g); MAP/SAS present under both labellings\n',
              suffix, nrow(chk), max(abs(chk$overall.x - chk$overall.y))))
}
scores = select(scores, -starts_with('n_metrics'))

# Kendall's W for m rankings of n items (columns = rankings), with the tie correction
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
kendall_w = function(r) {
  r = apply(r, 2, rank); m = ncol(r); n = nrow(r)
  S = sum((rowSums(r) - m * (n + 1) / 2)^2)
  Tj = sum(apply(r, 2, function(x) { t = table(x); sum(t^3 - t) }))
  12 * S / (m^2 * (n^3 - n) - m * Tj)
}

score_types = c(overall = 'Overall', bio = 'Bio conservation', batch = 'Batch correction')
long = scores |> pivot_longer(c(overall, bio, batch), names_to = 'score', values_to = 'value') |>
  pivot_wider(names_from = labelling, values_from = value)
paired = long |> filter(!is.na(best_annotator))
stopifnot(all(count(paired, tissue, case, score)$n == 23))

robust = paired |> group_by(tissue, case, score) |>
  summarise(n_methods = n(), tau = cor(ebcc, best_annotator, method = 'kendall'),
            weighted_tau = wtau(ebcc, best_annotator),   # higher score = more important
            W = kendall_w(cbind(-ebcc, -best_annotator)),
            top1_same = method[which.max(ebcc)] == method[which.max(best_annotator)],
            top3_overlap = length(intersect(method[order(-ebcc)][1:3], method[order(-best_annotator)][1:3])),
            .groups = 'drop')
write_csv(robust, file.path(out, paste0('label_robustness', suffix, '.csv')))

# consensus: mean rank over the 18 datasets under each labelling
ranks = paired |> group_by(tissue, case, score) |>
  mutate(rank_ebcc = rank(-ebcc), rank_alt = rank(-best_annotator)) |> ungroup()
mean_ranks = ranks |> group_by(score, method) |>
  summarise(rank_ebcc = mean(rank_ebcc), rank_alt = mean(rank_alt), .groups = 'drop')
consensus = mean_ranks |> group_by(score) |>
  summarise(tau = cor(rank_ebcc, rank_alt, method = 'kendall'), weighted_tau = wtau(-rank_ebcc, -rank_alt),
            W = kendall_w(cbind(rank_ebcc, rank_alt)), .groups = 'drop')
write_csv(mean_ranks, file.path(out, paste0('label_robustness_mean_ranks', suffix, '.csv')))

summ = robust |> group_by(score) |>
  summarise(tau_med = median(tau), tau_min = min(tau), tau_max = max(tau),
            wtau_med = median(weighted_tau), wtau_min = min(weighted_tau), wtau_max = max(weighted_tau),
            W_med = median(W), W_min = min(W), W_max = max(W),
            top1_same = sum(top1_same), n = n(), .groups = 'drop') |>
  left_join(consensus |> rename(tau_cons = tau, wtau_cons = weighted_tau, W_cons = W), by = 'score')
print(as.data.frame(summ), digits = 3)

# ---- figure ----
input_cols   = c('PCA / own' = '#2a78d6', 'UINMF' = '#eb6834', 'scVI latent' = '#1baf7a')
input_shapes = c('PCA / own' = 16, 'UINMF' = 17, 'scVI latent' = 15)
theme_bench = function(base_size = 11) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), plot.title.position = 'plot',
        strip.background = element_rect(fill = 'grey95', colour = 'grey70'),
        panel.border = element_rect(colour = 'grey70'))

# panel A: overall score only (the robustness question is about the overall ranking); one line per tissue
tissue_cols = c(bone_marrow = '#E69F00', fetal_liver = '#56B4E9', breast_cancer = '#CC79A7')   # as analysis/mcb_analysis.R
dot_d = robust |> filter(score == 'overall') |>
  mutate(tissue = factor(tissue, levels = names(tissue_cols)), case = factor(case, levels = case_levels),
         top1 = factor(ifelse(top1_same, 'Same', 'Different'), levels = c('Same', 'Different')))
ref = if (extra && file.exists(ref_file)) read_csv(ref_file, show_col_types = FALSE) else NULL
if (!is.null(ref) && !'statistic' %in% names(ref)) { message('reference file has no statistic column: rerun analysis/tau_reference.R, then this script'); ref = NULL }
if (!is.null(ref)) ref = ref |>
  filter(comparison == 'between datasets (EBCC)', score == 'overall', statistic == fig_stat)
stopifnot(is.null(ref) || nrow(ref) == 1)
dot_d = mutate(dot_d, stat = .data[[fig_stat]])
pA = ggplot(dot_d, aes(case, stat, colour = tissue, group = tissue)) +
  geom_blank() +                                   # discrete x before the reference layers
  { if (!is.null(ref)) list(
      annotate('rect', xmin = -Inf, xmax = Inf, ymin = ref$q25, ymax = ref$q75, fill = 'grey88'),
      geom_hline(yintercept = ref$median, linetype = 'dashed', colour = 'grey40', linewidth = 0.5),
      annotate('text', x = 5.45, y = ref$median, hjust = 1, vjust = -0.6, size = 3, colour = 'grey30',
               label = 'Between datasets of a tissue\n(median, interquartile range)')) } +
  geom_line(linewidth = 0.7, alpha = 0.8, position = position_dodge(width = 0.3)) +
  geom_point(aes(shape = top1), size = 3, stroke = 1, fill = 'white', position = position_dodge(width = 0.3)) +
  scale_colour_manual(values = tissue_cols, labels = tissue_labels[names(tissue_cols)], name = 'Tissue') +
  scale_shape_manual(values = c(Same = 16, Different = 21), name = 'Top-ranked method', drop = FALSE) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(add = 0.02)) +
  labs(x = 'Dataset', y = if (fig_stat == 'weighted_tau') expression("Weighted Kendall's" ~ tau[w] ~ '(EBCC vs best-annotator labels)') else
         expression("Kendall's" ~ tau ~ '(EBCC vs best-annotator labels)'),
       title = 'A  Agreement of the overall rankings, per dataset') +
  theme_bench() + theme(panel.grid.major.x = element_blank(), legend.position = 'bottom', legend.box = 'vertical',
                        legend.margin = margin(0, 0, 0, 0), axis.text.x = element_text(angle = 30, hjust = 1))

sc_d = mean_ranks |> filter(score == 'overall') |>
  mutate(input = case_when(str_detect(method, 'UINMF') ~ 'UINMF',
                           str_detect(method, 'scVI-latent') ~ 'scVI latent', TRUE ~ 'PCA / own'),
         input = factor(input, levels = names(input_cols)),
         base  = str_remove(method, ' (UINMF|scVI-latent)$'))
# colour = base method (legend ordered by its best mean rank), shape = input
base_order = sc_d |> group_by(base) |> summarise(r = min(rank_ebcc)) |> arrange(r) |> pull(base)
base_cols = setNames(c('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300',
                       '#6250d6', '#e34948', '#52514e', '#8a5a2b')[seq_along(base_order)], base_order)
sc_d = sc_d |> mutate(base = factor(base, levels = base_order))
cons_o = consensus |> filter(score == 'overall')
pB = ggplot(sc_d, aes(rank_ebcc, rank_alt)) +
  geom_abline(slope = 1, intercept = 0, colour = 'grey60', linetype = 'dashed') +
  geom_point(aes(fill = base, shape = input), size = 3.4, colour = 'white', stroke = 0.6) +
  annotate('text', x = 1, y = 23, hjust = 0, vjust = 1, size = 3.3,
           label = sprintf('%s == %.2f', if (fig_stat == 'weighted_tau') 'tau[w]' else 'tau', cons_o[[fig_stat]]), parse = TRUE) +
  scale_fill_manual(values = base_cols, name = 'Method',
                    guide = guide_legend(override.aes = list(shape = 21, size = 3.4), order = 1)) +
  scale_shape_manual(values = c('PCA / own' = 21, 'UINMF' = 24, 'scVI latent' = 22), name = 'Input',
                     guide = guide_legend(override.aes = list(fill = '#52514e', size = 3.4), order = 2)) +
  scale_x_continuous(limits = c(1, 23), breaks = c(1, 5, 10, 15, 20, 23)) +
  scale_y_continuous(limits = c(1, 23), breaks = c(1, 5, 10, 15, 20, 23)) +
  coord_equal() +
  labs(x = 'Mean rank, EBCC labels (1 = best)', y = 'Mean rank, best-annotator labels',
       title = sprintf('B  Mean overall rank over the %d datasets', n_distinct(paste(paired$tissue, paired$case)))) +
  theme_bench()

so = summ |> filter(score == 'overall')
w = fig_stat == 'weighted_tau'
st = if (w) "Weighted Kendall's $\\tau_w$" else "Kendall's $\\tau$"; sym = if (w) '$\\tau_w$' else '$\\tau$'
caption = sprintf(paste(
  'Robustness of the overall ranking to the cell-type labels (+MAP +SAS score, %d datasets: dissected and genes2--8 of bone marrow,',
  'fetal liver and breast cancer). (A) %s between the rankings of the 23 integration methods by overall score computed',
  'with EBCC labels and with the best single annotator\'s labels (per species, highest mean marker AUROC), for each dataset; open',
  'points: the top-ranked method differs between the two labellings. Median %s = %.2f (range %.3f--%.3f); same top-ranked method',
  'in %d of %d datasets. Grey band and dashed line: interquartile range and median of %s between the overall rankings',
  'of datasets of the same tissue under EBCC labels (%d pairs; median %.2f). (B) Mean rank of each method over the %d datasets',
  'under the two labellings (%s = %.2f); dashed line: equal rank. Colon is not included: EBCC is the best annotator in both',
  'species there, so the two labellings coincide. scANVI is trained on EBCC labels under both labellings.'),
  so$n, st, sym, if (w) so$wtau_med else so$tau_med, if (w) so$wtau_min else so$tau_min, if (w) so$wtau_max else so$tau_max,
  so$top1_same, so$n, sym, ref$n, ref$median, so$n, sym, if (w) so$wtau_cons else so$tau_cons)
writeLines(caption, file.path(out, paste0('12_label_robustness', suffix, '_caption.txt')))

p = pA + pB + plot_layout(widths = c(1, 1.5))
ggsave(file.path(out, paste0('12_label_robustness', suffix, '.png')), p, width = 15, height = 7.5, dpi = 300, bg = 'white')
ggsave(file.path(out, paste0('12_label_robustness', suffix, '.pdf')), p, width = 15, height = 7.5, bg = 'white', device = cairo_pdf)
ggsave(file.path(out, paste0('12a_label_robustness_overall', suffix, '.png')), pA, width = 5.5, height = 6, dpi = 300, bg = 'white')
ggsave(file.path(out, paste0('12b_label_robustness_scatter', suffix, '.png')), pB, width = 8.5, height = 7, dpi = 300, bg = 'white')
