# Main benchmark figure, appendix heatmap and orthologue-removal figure:  Rscript analysis/benchmark_figures.R
# Scores as in analysis/benchmark_summary.R (unchanged): batch = mean of PCR_batch, ASW_label/batch, iLISI, graph_conn;
# bio = mean of NMI, ARI, ASW_label, cLISI; overall = 0.6 bio + 0.4 batch. Graph-only methods (BBKNN x3, SAMap)
# have no embedding, so PCR_batch, ASW_label/batch and ASW_label are missing: their batch is the mean of
# iLISI and graph_conn, their bio the mean of NMI, ARI and cLISI.
# Improvement = method score - score of 'Unintegrated' (joint log-normalised PCA, 05.00_unintegrated.csv) in the
# same tissue and dataset. Bootstrap SD of the improvement: from the paired deltas in bootstrap/summary/1_long.csv
# (50 subsamples of 20,000 cells, dissected); not valid for graph-only methods (their stored graph cannot be
# subsampled), so left empty for them.
# Variant (config extra_bio / extra_batch): the extra metrics from analysis/extra_metrics.py (MAP, NC) and analysis/sas_scores.py (SAS)
# are added as further equally weighted metrics of the bio / batch means; the scores above stay the reference
# definition. Graph-only methods have no embedding, so no MAP / NC / SAS: their scores are unchanged.
# Output (figures/benchmark_summary/): benchmark_improvement<suffix>.{pdf,png,csv}, benchmark_strategy<suffix>.{pdf,png},
#   benchmark_heatmap<suffix>.{pdf,png},
#   orthologue_removal<suffix>.{pdf,png}; printout on stdout. suffix = '' without extra metrics, else e.g. '_MAP_SAS'.
source('library.R')

# ---------------- config ----------------
out          = 'figures/benchmark_summary'
boot_file    = 'bootstrap/summary/1_long.csv'
reference    = 'Unintegrated'                                   # native input, log-normalised joint PCA
ref_rows     = c('Unintegrated', 'Unintegrated scVI-latent', 'Unintegrated UINMF')   # heatmap reference rows
graph_only   = c('BBKNN', 'BBKNN UINMF', 'BBKNN scVI-latent', 'SAMap')
label_based  = 'scANVI'                                         # trained on EBCC labels: marked with a dagger
cases        = c(dissected = 0, genes2 = 20, genes4 = 40, genes6 = 60, genes8 = 80)   # % orthologues removed
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
tissue_cols  = c(bone_marrow = '#E69F00', colon = '#009E73', breast_cancer = '#CC79A7', fetal_liver = '#56B4E9')  # Okabe-Ito, as analysis/ebcc_support.R
tissue_shapes = c(bone_marrow = 16, colon = 17, breast_cancer = 15, fetal_liver = 18)
input_levels = c('Native', 'scVI-latent', 'UINMF')
input_fill   = c(Native = '#3a3a3a', `scVI-latent` = '#8a8a8a', UINMF = '#d0d0d0')
input_text   = c(Native = 'N', `scVI-latent` = 'S', UINMF = 'U')
# integration strategy of the algorithm (a method on UINMF / scVI-latent input keeps its algorithm's strategy)
strategy     = c(FastMNN = 'match', Scanorama = 'match', CCA = 'match', RPCA = 'match', SAMap = 'match',
                 BBKNN = 'balance', Harmony = 'balance',
                 scVI = 'joint', scANVI = 'joint', LIGER = 'joint')
strategy_labels = c(match = 'Matching corresponding cells', balance = 'Balancing species in a shared space',
                    joint = 'Learning a joint representation')
input_hl     = c(Native = '#d9d9d9', `scVI-latent` = '#7B5EA7', UINMF = '#CCBB44')   # improvement figure: row tint (none for native), clear of the tissue colours
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI')
extra_bio     = c('MAP')        # added to bio: any of 'MAP', 'NC'; character(0) for the reference definition
extra_batch   = c('SAS')        # added to batch: 'SAS' or character(0)
compare_extra = list(`+MAP +SAS` = list(bio = 'MAP', batch = 'SAS'), `+MAP +NC +SAS` = list(bio = c('MAP', 'NC'), batch = 'SAS'),
                     `+MAP` = list(bio = 'MAP', batch = character(0)), `+SAS` = list(bio = character(0), batch = 'SAS'))  # printout only
extra_file    = 'extra_metrics/extra_metrics_all.csv'
sas_file      = 'sas/sas_all.csv'
boot_extra_file = 'bootstrap_extra/extra_long.csv'
graph_metrics_only = FALSE      # TRUE: every method scored only on the metrics graph-only methods also have (bio: NMI, ARI,
                                # cLISI; batch: iLISI, graph_conn), no extra metrics; same weights; suffix '_graph_metrics'
# -----------------------------------------
bio_metrics_ref = bio_metrics; batch_metrics_ref = batch_metrics   # reference definition, for the bootstrap-table check
if (graph_metrics_only) {
  bio_metrics = c('NMI_cluster/label', 'ARI_cluster/label', 'cLISI'); batch_metrics = c('iLISI', 'graph_conn')
  extra_bio = character(0); extra_batch = character(0)
}

score_dir = function(tissue, case) {   # as analysis/benchmark_summary.R
  if (tissue == 'breast_cancer')
    return(sprintf('breast_cancer/data/%s/', if (case %in% c('full', 'dissected')) paste0(case, '_mouse') else case))
  if (case %in% c('full', 'dissected')) sprintf('%s/data/05_evaluated_score/%s_mouse/', tissue, case)
  else sprintf('%s/data/%s/', tissue, case)
}
read_metrics = function(f) {
  df = suppressMessages(read_csv(f, show_col_types = FALSE))
  tibble(metric = df[[1]], value = as.numeric(df[[2]]), method = names(df)[2])
}
input_of = function(m) factor(case_when(str_detect(m, 'UINMF') ~ 'UINMF', str_detect(m, 'scVI-latent') ~ 'scVI-latent',
                                        TRUE ~ 'Native'), levels = input_levels)
label_of = function(m) case_when(m %in% c('CCA', 'RPCA') ~ paste0(m, ' (Seurat)'),
                                 m %in% c('CCA UINMF', 'RPCA UINMF', 'CCA scVI-latent', 'RPCA scVI-latent') ~
                                   str_replace(m, '^(CCA|RPCA)', '\\1 (Seurat)'),
                                 m == label_based ~ paste0(m, ' †'),
                                 m %in% graph_only ~ paste0(m, ' ‡'), TRUE ~ m)

# ---- full-data scores, all cases ----
suffix = if (graph_metrics_only) '_graph_metrics' else if (length(c(extra_bio, extra_batch))) paste0('_', paste(c(extra_bio, extra_batch), collapse = '_')) else ''
met = expand_grid(tissue = names(tissue_labels), case = names(cases)) |>
  pmap_dfr(function(tissue, case)
    list.files(score_dir(tissue, case), pattern = '^05\\..*\\.csv$', full.names = TRUE) |>
      map_dfr(read_metrics) |> mutate(tissue = tissue, case = case)) |>
  filter(!str_detect(method, 'scIB'), !is.na(value))
extra = read_csv(extra_file, show_col_types = FALSE) |> select(tissue, case = dataset, method, MAP, NC) |>
  full_join(read_csv(sas_file, show_col_types = FALSE) |> select(tissue, case = dataset, method, SAS = sas),
            by = c('tissue', 'case', 'method')) |>
  mutate(case = str_remove(case, '_mouse$')) |> filter(case %in% names(cases)) |>
  pivot_longer(c(MAP, NC, SAS), names_to = 'metric', values_to = 'value') |> filter(!is.na(value))
met = bind_rows(met, extra)
# bio / batch = equal-weight means of the listed metrics that exist for the method; overall = 0.6 bio + 0.4 batch
score_with = function(m, bio_m, batch_m, bio_base = bio_metrics, batch_base = batch_metrics) m |>
  mutate(group = case_when(metric %in% c(batch_base, batch_m) ~ 'batch', metric %in% c(bio_base, bio_m) ~ 'bio')) |>
  filter(!is.na(group)) |>
  group_by(tissue, case, method, group) |> summarise(score = mean(value), n_metrics = n(), .groups = 'drop') |>
  pivot_wider(names_from = group, values_from = c(score, n_metrics)) |>
  rename(bio = score_bio, batch = score_batch) |>
  mutate(overall = 0.6 * bio + 0.4 * batch)
improve = function(sc) {
  ref = sc |> filter(method == reference) |> select(tissue, case, ref_overall = overall, ref_bio = bio, ref_batch = batch)
  stopifnot(nrow(ref) == length(cases) * 4)
  sc |> filter(!method %in% ref_rows) |> left_join(ref, by = c('tissue', 'case')) |>
    mutate(imp_overall = overall - ref_overall, imp_bio = bio - ref_bio, imp_batch = batch - ref_batch,
           input = input_of(method))
}
scores0 = score_with(met, character(0), character(0))           # reference definition
scores  = score_with(met, extra_bio, extra_batch)               # this run
imp0 = improve(scores0); imp = improve(scores)
# MAP / NC / SAS must exist for every embedding method and reference in the datasets used
miss = met |> filter(metric %in% c(extra_bio, extra_batch)) |> distinct(tissue, case, method, metric) |>
  right_join(expand_grid(distinct(filter(scores, !method %in% graph_only), tissue, case, method),
                         metric = c(extra_bio, extra_batch)), by = c('tissue', 'case', 'method', 'metric')) |>
  anti_join(filter(met, metric %in% c(extra_bio, extra_batch)), by = c('tissue', 'case', 'method', 'metric'))
missing_cases = unique(miss$case)
if (length(missing_cases)) cat('extra metrics missing for', nrow(miss), 'method x dataset x metric rows, cases:',
                               paste(missing_cases, collapse = ', '), '-> those cases are skipped for this variant\n')
stopifnot(!'dissected' %in% missing_cases)

# ---- bootstrap SD of the paired improvement (dissected) ----
boot = read_csv(boot_file, show_col_types = FALSE)
# check: the bootstrap table's full-data scores and deltas reproduce the score files (reference definition)
chk = boot |> filter(`repeat` == 'point', metric %in% c('overall', 'bio', 'batch')) |>
  select(tissue = dataset, method, metric, score, delta = delta_from_unintegrated) |>
  inner_join(improve(score_with(met, character(0), character(0), bio_metrics_ref, batch_metrics_ref)) |> filter(case == 'dissected') |>
               pivot_longer(c(overall, bio, batch), names_to = 'metric', values_to = 'file_score') |>
               mutate(file_delta = file_score - case_when(metric == 'overall' ~ ref_overall, metric == 'bio' ~ ref_bio,
                                                         TRUE ~ ref_batch)) |>
               select(tissue, method, metric, file_score, file_delta), by = c('tissue', 'method', 'metric'))
cat(sprintf('check vs bootstrap table: %d rows, max |score diff| %.2g, max |delta diff| %.2g\n',
            nrow(chk), max(abs(chk$score - chk$file_score)), max(abs(chk$delta - chk$file_delta))))
stopifnot(nrow(chk) > 0, max(abs(chk$score - chk$file_score)) < 1e-6, max(abs(chk$delta - chk$file_delta)) < 1e-6)
# per-subsample metrics of every method and of the reference (reference value = score - paired delta)
bm = boot |> filter(`repeat` != 'point', metric %in% c(batch_metrics, bio_metrics), !method %in% graph_only) |>
  transmute(tissue = dataset, rep = as.integer(`repeat`), method, metric, value = score,
            ref_value = score - delta_from_unintegrated)
bm_ref = bm |> group_by(tissue, rep, metric) |> summarise(spread = max(ref_value) - min(ref_value), value = mean(ref_value),
                                                           .groups = 'drop') |> mutate(method = reference)
stopifnot(max(bm_ref$spread) < 1e-9)   # the same reference value for every method of a subsample
bm_ref = select(bm_ref, -spread)
bm = bind_rows(select(bm, -ref_value), bm_ref)
if (length(c(extra_bio, extra_batch)))
  bm = bind_rows(bm, read_csv(boot_extra_file, show_col_types = FALSE) |>
                   pivot_longer(c(MAP, NC, SAS), names_to = 'metric', values_to = 'value') |>
                   filter(metric %in% c(extra_bio, extra_batch), method %in% unique(bm$method)) |>
                   select(tissue, rep, method, metric, value))
bs = bm |> rename(case = rep) |> score_with(extra_bio, extra_batch)
bref = bs |> filter(method == reference) |> select(tissue, case, r_o = overall, r_bi = bio, r_ba = batch)
sds = bs |> filter(method != reference) |> inner_join(bref, by = c('tissue', 'case')) |>
  group_by(tissue, method) |>
  summarise(boot_sd_overall = sd(overall - r_o), boot_sd_bio = sd(bio - r_bi), boot_sd_batch = sd(batch - r_ba),
            boot_n = n(), .groups = 'drop')

dis = imp |> filter(case == 'dissected') |> left_join(sds, by = c('tissue', 'method'))
write_csv(dis |> transmute(method, input, graph_only = method %in% graph_only, tissue, overall, bio, batch,
                           ref_overall, ref_bio, ref_batch, imp_overall, imp_bio, imp_batch,
                           boot_sd_overall, boot_sd_bio, boot_sd_batch, boot_n, n_metrics_bio, n_metrics_batch) |>
            arrange(method, tissue), file.path(out, paste0('benchmark_improvement', suffix, '.csv')))

# ---- printout ----
summ = dis |> group_by(method, input) |>
  summarise(mean_imp_overall = mean(imp_overall), mean_imp_bio = mean(imp_bio), mean_imp_batch = mean(imp_batch),
            beats_overall = sum(imp_overall > 0), beats_bio = sum(imp_bio > 0), beats_batch = sum(imp_batch > 0),
            .groups = 'drop')
rk = dis |> group_by(tissue) |> mutate(r = rank(-overall)) |> group_by(method) |> summarise(mean_rank = mean(r))
summ = summ |> left_join(rk, by = 'method') |> arrange(-mean_imp_overall) |>
  mutate(pos_score = row_number(), pos_rank = rank(mean_rank, ties.method = 'min'),
         flag = ifelse(abs(pos_score - pos_rank) > 2, '<-- moves >2 places', ''))
cat('\n== Dissected: mean improvement over Unintegrated (native) and tissues beating it (of 4) ==\n')
summ |> transmute(method, input, overall = round(mean_imp_overall, 3), bio = round(mean_imp_bio, 3),
                  batch = round(mean_imp_batch, 3), `beats o/b/b` = paste(beats_overall, beats_bio, beats_batch, sep = '/'),
                  pos_score, mean_rank = round(mean_rank, 2), pos_rank, flag) |>
  as.data.frame() |> print(row.names = FALSE)
cat('\n== Bootstrap SD of the improvement (embedding methods, per method x tissue) ==\n')
dis |> filter(!method %in% graph_only) |>
  pivot_longer(c(boot_sd_overall, boot_sd_bio, boot_sd_batch), names_to = 'score', values_to = 'sd') |>
  group_by(score = str_remove(score, 'boot_sd_')) |>
  summarise(median = median(sd), min = min(sd), max = max(sd), n = n()) |>
  mutate(across(c(median, min, max), \(x) round(x, 4))) |> as.data.frame() |> print(row.names = FALSE)
cat('\nBy tissue (overall):\n')
dis |> filter(!method %in% graph_only) |> group_by(tissue) |>
  summarise(median = round(median(boot_sd_overall), 4), min = round(min(boot_sd_overall), 4),
            max = round(max(boot_sd_overall), 4)) |> as.data.frame() |> print(row.names = FALSE)

# ranking under the reference definition and with the extra metrics (dissected, mean overall improvement)
positions = function(sc, name) improve(sc) |> filter(case == 'dissected') |> group_by(method) |>
  summarise(v = mean(imp_overall)) |> mutate(variant = name, pos = rank(-v, ties.method = 'min'))
pos = bind_rows(positions(scores0, 'reference'),
                imap_dfr(compare_extra, \(e, nm) positions(score_with(met, e$bio, e$batch), nm)))
cat('\n== Position by mean overall improvement (dissected): reference definition vs with extra metrics ==\n')
pw = pos |> select(method, variant, pos) |> pivot_wider(names_from = variant, values_from = pos) |> arrange(reference)
print(as.data.frame(pw), row.names = FALSE)
emb = filter(pos, !method %in% graph_only)
for (nm in names(compare_extra)) {
  x = inner_join(filter(emb, variant == 'reference'), filter(emb, variant == nm), by = 'method')
  cat(sprintf('%-15s Kendall tau vs reference (embedding methods) %.2f; top 5 %s; moved >2 places: %s\n', nm,
              cor(x$v.x, x$v.y, method = 'kendall'),
              paste(x$method[order(x$pos.y)][1:5], collapse = ', '),
              paste(x$method[abs(x$pos.x - x$pos.y) > 2], collapse = ', ')))
}
cat('\nGraph-only methods: scores over the available metrics -- batch from',
    paste(unique(dis$n_metrics_batch[dis$method %in% graph_only]), collapse = '/'), 'metrics (iLISI, graph_conn), bio from',
    paste(unique(dis$n_metrics_bio[dis$method %in% graph_only]), collapse = '/'), 'metrics (NMI, ARI, cLISI); no bootstrap SD.\n')

# ---- 1. main figure: improvement over unintegrated ----
ord = summ$method   # best first
y_levels = rev(label_of(ord))
long = dis |> pivot_longer(c(imp_overall, imp_bio, imp_batch), names_to = 'score', values_to = 'imp') |>
  mutate(score = factor(str_remove(score, 'imp_'), levels = c('overall', 'bio', 'batch'),
                        labels = c('Overall (0.6 bio + 0.4 batch)', 'Bio-conservation', 'Batch correction')),
         y = factor(label_of(method), levels = y_levels),
         tissue = factor(tissue, levels = names(tissue_labels)))
long = long |> mutate(strategy = factor(strategy[word(method, 1)], levels = names(strategy_labels), labels = strategy_labels))
means = long |> group_by(method, y, strategy, score) |> summarise(imp = mean(imp), .groups = 'drop')
strip = tibble(method = ord, y = factor(label_of(ord), levels = y_levels), input = input_of(ord),
               strategy = factor(strategy[word(ord, 1)], levels = names(strategy_labels), labels = strategy_labels))
stopifnot(!is.na(strip$strategy), !is.na(long$strategy))
position = setNames(seq_along(ord), label_of(ord))   # place by mean overall improvement, shown next to the name

panel = function(s) {
  ggplot(filter(long, score == s), aes(imp, y)) +
    geom_blank() +   # sets the row order before the shading layer
    geom_tile(data = filter(strip, input != 'Native'), aes(x = 0, y = y, fill = input), width = Inf, height = 0.92,
              alpha = 0.22, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey40', linewidth = 0.4) +
    geom_point(aes(colour = tissue, shape = tissue), size = 2.3, alpha = 0.9) +
    geom_point(data = filter(means, score == s), shape = '|', size = 5.5, colour = 'black') +
    scale_colour_manual(values = tissue_cols, labels = tissue_labels, name = 'Tissue') +
    scale_shape_manual(values = tissue_shapes, labels = tissue_labels, name = 'Tissue') +
    scale_fill_manual(values = input_hl, guide = 'none') +
    labs(x = 'Improvement over unintegrated', y = NULL, title = s) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.title = element_text(face = 'bold', size = 11, hjust = 0.5),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
strip_p = ggplot(strip, aes(1, y)) +
  geom_text(aes(x = 0, label = position[as.character(y)]), size = 3.4, fontface = 'bold', colour = 'grey10') +
  geom_tile(aes(fill = input), width = 0.9, height = 0.85) +
  geom_text(aes(label = input_text[input], colour = input == 'Native'), size = 2.6, fontface = 'bold') +
  scale_x_continuous(limits = c(-0.5, 1.5)) +
  scale_fill_manual(values = input_hl, name = 'Input', drop = FALSE, labels = \(x) sprintf('%s (%s)', x, input_text[x])) +
  scale_colour_manual(values = c(`TRUE` = 'grey20', `FALSE` = 'white'), guide = 'none') +
  labs(title = ' ') +
  theme_void(base_size = 11) +
  theme(axis.text.y = element_text(colour = 'grey10', hjust = 1, size = 9.5, margin = margin(r = 2)),
        plot.title = element_text(size = 11))
metric_names = c(`NMI_cluster/label` = 'NMI', `ARI_cluster/label` = 'ARI', ASW_label = 'ASW label', cLISI = 'cLISI',
                 PCR_batch = 'PCR batch', `ASW_label/batch` = 'ASW batch', iLISI = 'iLISI', graph_conn = 'graph connectivity',
                 MAP = 'MAP', NC = 'NC', SAS = 'SAS')
mlist = function(x) paste(metric_names[x], collapse = ', ')
top_n = 7
gaps = summ$mean_imp_overall[1:top_n]
cap1 = paste0(
  'Point: one tissue; black bar: mean over the four tissues. Improvement = method score - score of unintegrated data ',
  '(log-normalised joint PCA) in the same tissue.\nOverall = 0.6 bio-conservation + 0.4 batch correction. ',
  'Bio-conservation = mean of ', mlist(c(bio_metrics, extra_bio)), '.\nBatch correction = mean of ',
  mlist(c(batch_metrics, extra_batch)), '.\n',
  'Number: place by mean overall improvement. Ranks 2-', top_n, ' lie within ', sprintf('%.3f', gaps[2] - gaps[top_n]),
  ';\nthe median bootstrap SD of a method\'s improvement (50 subsamples of 20,000 cells) is ',
  sprintf('%.3f', median(dis$boot_sd_overall, na.rm = TRUE)),
  if (graph_metrics_only) ' (embedding methods).\n' else ', so neighbouring places are not distinguishable.\n',
  'Input: scVI and scANVI start from counts, so they are native; their latent space is the scVI-latent input of other methods. ',
  '\u2020 trained on the EBCC labels.\n',
  if (graph_metrics_only) paste0('\u2021 graph-based (', paste(graph_only, collapse = ', '), '): no embedding, so no bootstrap SD. Every method here is scored only on the ',
                                 'metrics these methods also have,\nso all rows use the same metrics.')
  else paste0('\u2021 graph-based (', paste(graph_only, collapse = ', '), '): no embedding, so only ', mlist(intersect(bio_metrics, c('NMI_cluster/label', 'ARI_cluster/label', 'cLISI'))),
  ' (bio) and ', mlist(c('iLISI', 'graph_conn')), ' (batch) apply; their scores are the means of these.'))
p1 = strip_p + panel(levels(long$score)[1]) + panel(levels(long$score)[2]) + panel(levels(long$score)[3]) +
  patchwork::plot_layout(widths = c(0.16, 1.6, 1, 1), guides = 'collect') +
  patchwork::plot_annotation(caption = cap1, theme = theme(plot.caption = element_text(hjust = 0, colour = 'grey30', size = 9.5,
                                                                                         lineheight = 1.15))) &
  theme(legend.position = 'bottom')
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('benchmark_improvement', suffix, '.', ext)), p1, width = 13, height = 9, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')

# ---- 1a. improvement per metric (dissected), rows and order as figure 1 ----
# each metric of this run's bio / batch definition: method value - Unintegrated value in the same tissue
per_metric = c(bio_metrics, extra_bio, batch_metrics, extra_batch)
pm = met |> filter(case == 'dissected', metric %in% per_metric)
pm_ref = pm |> filter(method == reference) |> select(tissue, metric, ref = value)
pm = pm |> filter(method %in% ord) |> inner_join(pm_ref, by = c('tissue', 'metric')) |>
  mutate(imp = value - ref, y = factor(label_of(method), levels = y_levels), tissue = factor(tissue, levels = names(tissue_labels)),
         group = ifelse(metric %in% c(bio_metrics, extra_bio), 'Bio', 'Batch'),
         metric = factor(metric, levels = per_metric, labels = paste0(ifelse(per_metric %in% c(bio_metrics, extra_bio), 'Bio: ', 'Batch: '),
                                                                     metric_names[per_metric])))
write_csv(pm |> select(method, tissue, metric, value, ref, imp) |> arrange(metric, method, tissue),
          file.path(out, paste0('benchmark_metrics', suffix, '.csv')))
pm_means = pm |> group_by(metric, y) |> summarise(imp = mean(imp), .groups = 'drop')
pm_panel = function(m) {
  d = filter(pm, metric == m)
  ggplot(d, aes(imp, y)) +
    geom_blank(data = strip, aes(x = 0, y = y), inherit.aes = FALSE) +
    geom_tile(data = filter(strip, input != 'Native'), aes(x = 0, y = y, fill = input), width = Inf, height = 0.92,
              alpha = 0.22, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey40', linewidth = 0.4) +
    geom_point(aes(colour = tissue, shape = tissue), size = 2.1, alpha = 0.9) +
    geom_point(data = filter(pm_means, metric == m), shape = '|', size = 5.5, colour = 'black') +
    scale_colour_manual(values = tissue_cols, labels = tissue_labels, name = 'Tissue') +
    scale_shape_manual(values = tissue_shapes, labels = tissue_labels, name = 'Tissue') +
    scale_fill_manual(values = input_hl, guide = 'none') +
    scale_y_discrete(drop = FALSE) +
    labs(x = 'Improvement over unintegrated', y = NULL, title = m) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.title = element_text(face = 'bold', size = 11, hjust = 0.5),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(5.5, 14, 5.5, 5.5))
}
# one row of panels for bio metrics, one for batch, each with the method strip; empty columns pad the shorter row
mrow = function(ms, n_col) {
  ps = c(list(strip_p), lapply(ms, pm_panel), rep(list(patchwork::plot_spacer()), n_col - length(ms)))
  patchwork::wrap_plots(ps, nrow = 1, widths = c(0.16 * 4 / n_col * 1.6, rep(1, n_col)))
}
lv = levels(pm$metric); bio_lv = lv[startsWith(lv, 'Bio')]; batch_lv = lv[startsWith(lv, 'Batch')]
n_col = max(length(bio_lv), length(batch_lv))
p1a = (mrow(bio_lv, n_col) / mrow(batch_lv, n_col)) + patchwork::plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('benchmark_metrics', suffix, '.', ext)), p1a, width = 3 + 2.6 * n_col, height = 15.5,
         dpi = 300, bg = 'white', device = if (ext == 'pdf') cairo_pdf else 'png')
writeLines(paste0('Improvement of each metric over unintegrated data (log-normalised joint PCA) in the same tissue (dissected mouse data). ',
                  'Point: one tissue; black bar: mean over the four tissues. Rows, order and numbers as in benchmark_improvement', suffix,
                  ' (place by mean overall improvement). Shading: input (scVI-latent, UINMF). \u2020 trained on the EBCC labels; ',
                  '\u2021 graph-based (no embedding).',
                  if (!graph_metrics_only) ' Metrics a graph-based method lacks are left empty in its row.' else ''),
           file.path(out, paste0('benchmark_metrics', suffix, '_caption.txt')))

# ---- 1b. which strategy works best, on which input ----
# Rows = strategy, columns = input; small points = the algorithms of that cell (mean overall improvement over the
# four tissues, as the bar in figure 1); diamond and number = their mean, n = algorithms in the cell.
st = summ |> transmute(method, input, imp = mean_imp_overall, algorithm = word(method, 1),
                       strategy = factor(strategy[algorithm], levels = rev(names(strategy_labels)), labels = rev(strategy_labels)),
                       algo_lab = case_when(algorithm %in% c('CCA', 'RPCA') ~ paste0(algorithm, ' (Seurat)'), TRUE ~ algorithm),
                       algo_lab = paste0(algo_lab, ifelse(method %in% graph_only, ' ‡', ''), ifelse(method == label_based, ' †', '')))
stopifnot(!is.na(st$strategy))
# strategy mean only over algorithms run on every input of their strategy (paired), so the inputs compare like with like;
# no mean when fewer than 2 such algorithms (joint representation: LIGER only)
st = st |> group_by(strategy) |> mutate(n_inputs = n_distinct(input)) |> group_by(strategy, algorithm) |>
  mutate(paired = n_distinct(input) == n_inputs) |> ungroup()
st_mean = st |> filter(paired) |> group_by(strategy) |> filter(n_distinct(algorithm) >= 2) |>
  group_by(strategy, input) |> summarise(imp = mean(imp), n = n(), .groups = 'drop')
cat('\n== Mean overall improvement by strategy x input (paired algorithms only) ==\n')
st_mean |> mutate(txt = sprintf('%.3f (n=%d)', imp, n)) |> select(strategy, input, txt) |>
  pivot_wider(names_from = input, values_from = txt) |> as.data.frame() |> print(row.names = FALSE)
write_csv(st |> select(strategy, input, method, imp), file.path(out, paste0('benchmark_strategy', suffix, '.csv')))
# rows = algorithms grouped by strategy (best first), plus a shaded 'Strategy mean' row on top of each group
mean_lab = st_mean |> distinct(strategy, n) |> transmute(strategy, lab = sprintf('Paired mean (n = %d)', n))
st_mean = st_mean |> left_join(mean_lab, by = 'strategy')
st_rows = bind_rows(st |> select(strategy, input, row = algo_lab, imp),
                    st_mean |> transmute(strategy, input, row = lab, imp)) |>
  mutate(is_mean = row %in% mean_lab$lab,
         offset = c(Native = 0.15, `scVI-latent` = 0, UINMF = -0.15)[as.character(input)])   # keeps equal values apart
# row position within each strategy panel: mean row on top, then algorithms by their best value
row_pos = st |> group_by(strategy, row = algo_lab) |> summarise(best = max(imp), .groups = 'drop') |>
  group_by(strategy) |> arrange(best, .by_group = TRUE) |> mutate(yr = row_number()) |>
  bind_rows(x = _, mean_lab |> rename(row = lab)) |>
  group_by(strategy) |> mutate(yr = ifelse(row %in% mean_lab$lab, max(yr, na.rm = TRUE) + 1.2, yr)) |> ungroup() |> select(-best)
st_rows = st_rows |> left_join(row_pos, by = c('strategy', 'row')) |>
  mutate(y = yr + offset, strategy = factor(strategy, levels = rev(levels(strategy))))
row_pos = row_pos |> mutate(strategy = factor(strategy, levels = levels(st_rows$strategy)), is_mean = row %in% mean_lab$lab)
x_left = min(st_rows$imp) - 0.006
p1b = ggplot(st_rows, aes(imp, y)) +
  geom_rect(data = filter(row_pos, is_mean), aes(xmin = -Inf, xmax = Inf, ymin = yr - 0.5, ymax = yr + 0.5),
            fill = 'grey92', inherit.aes = FALSE) +
  geom_line(aes(y = yr, group = row), colour = 'grey65', linewidth = 0.5) +
  geom_point(data = filter(st_rows, !is_mean), aes(fill = input, shape = 'alg'), size = 3, colour = 'grey20', stroke = 0.4) +
  geom_point(data = filter(st_rows, is_mean), aes(fill = input, shape = 'mean'), size = 4.6, colour = 'black', stroke = 0.7) +
  geom_text(data = row_pos, aes(x = -Inf, y = yr, label = paste0(row, '  '), fontface = ifelse(is_mean, 'bold', 'plain')),
            hjust = 1, nudge_x = 0, size = 3.5, colour = 'grey10', inherit.aes = FALSE) +
  facet_grid(strategy ~ ., scales = 'free_y', space = 'free_y', labeller = labeller(strategy = label_wrap_gen(16))) +
  coord_cartesian(clip = 'off') +
  scale_fill_manual(values = input_hl, name = 'Input', guide = guide_legend(order = 1, override.aes = list(shape = 22, size = 5))) +
  scale_shape_manual(values = c(alg = 21, mean = 23), name = NULL,
                     labels = c(alg = 'One algorithm', mean = 'Paired mean of the strategy\n(algorithms run on every input)'),
                     guide = guide_legend(order = 2, override.aes = list(fill = 'grey60', size = c(3, 4.6)))) +
  scale_x_continuous(breaks = scales::breaks_pretty(6)) +
  scale_y_continuous(expand = expansion(add = 0.6)) +
  labs(x = 'Mean overall improvement over unintegrated (four tissues)', y = NULL,
       caption = paste0('SAMap, scVI and scANVI (one input) are shown but not averaged; ',
                        'joint representation has one paired algorithm (LIGER), so no mean.\n',
                        'A method on UINMF / scVI-latent input keeps its ',
                        'algorithm\'s strategy.\nscVI-latent is not used with joint-representation methods (scVI is its source). ',
                        '‡ graph-only, † trained on labels.')) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.title.position = 'plot',
        plot.title = element_text(face = 'bold'), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 125),
        strip.background = element_rect(fill = 'grey95', colour = 'grey70'),
        strip.text.y = element_text(face = 'bold', angle = 0, hjust = 0),
        legend.position = 'top', plot.caption = element_text(hjust = 0, colour = 'grey35'), plot.caption.position = 'plot')
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('benchmark_strategy', suffix, '.', ext)), p1b, width = 9.5, height = 7, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')

# ---- 1c. consistency of the ranking across tissues: pairwise weighted Kendall tau, Kendall's W ----
# Methods ranked within each tissue by overall improvement (= by overall score: the reference is fixed within a tissue),
# all methods of figure 1. W with a permutation p-value (ranks shuffled independently within each tissue).
# Figure: weighted Kendall tau (wtau below); the plain Kendall tau is still printed and saved in the csv.
n_perm = 10000
rk_t = dis |> select(tissue, method, imp_overall) |> group_by(tissue) |> mutate(r = rank(-imp_overall)) |> ungroup() |>
  select(-imp_overall) |> pivot_wider(names_from = tissue, values_from = r)
rmat = as.matrix(rk_t[, names(tissue_labels)])
kendall_w = function(m) { k = ncol(m); n = nrow(m); R = rowSums(m); 12 * sum((R - mean(R))^2) / (k^2 * (n^3 - n)) }
W = kendall_w(rmat)
set.seed(1)
W_perm = replicate(n_perm, kendall_w(apply(rmat, 2, sample)))
W_p = (1 + sum(W_perm >= W)) / (n_perm + 1)
# Weighted Kendall tau (Vigna 2015): hyperbolic weight 1/(r + 1) for the item at importance rank r (0 = best), a pair weighs
# the sum of its two weights, so swaps near the top count most; averaged over both rankings as the importance order.
# Equals scipy.stats.weightedtau(x, y) with its defaults (checked to 1e-15, with and without ties). x, y: higher = better.
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
tau = expand_grid(a = names(tissue_labels), b = names(tissue_labels)) |>
  mutate(tau = map2_dbl(a, b, \(x, y) cor(rmat[, x], rmat[, y], method = 'kendall')),
         wtau = map2_dbl(a, b, \(x, y) wtau(-rmat[, x], -rmat[, y])),     # rank 1 = best = most important
         ia = match(a, names(tissue_labels)), ib = match(b, names(tissue_labels)))
cat(sprintf('\n== Ranking consistency across tissues (%d methods, overall): Kendall W = %.2f, permutation p = %.4f ==\n',
            nrow(rmat), W, W_p))
tau |> filter(ia < ib) |> transmute(pair = paste(a, b, sep = ' vs '), tau = round(tau, 2), weighted_tau = round(wtau, 2)) |>
  as.data.frame() |> print(row.names = FALSE)
write_csv(tau |> filter(ia != ib) |> select(a, b, tau, weighted_tau = wtau), file.path(out, paste0('benchmark_tissue_tau', suffix, '.csv')))
tl2 = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast\ncancer', fetal_liver = 'Fetal liver')
tau_off = filter(tau, ia > ib)   # lower triangle; diagonal (tau = 1) shown in grey, so all four tissues get a row and a column
pC = ggplot(tau, aes(factor(tl2[b], levels = tl2), factor(tl2[a], levels = rev(tl2)))) +
  geom_tile(data = filter(tau, ia == ib), fill = 'grey85', colour = 'white', linewidth = 1) +
  geom_text(data = filter(tau, ia == ib), label = '1', colour = 'grey45', size = 4.2) +
  geom_tile(data = tau_off, aes(fill = wtau), colour = 'white', linewidth = 1) +
  geom_text(data = tau_off, aes(label = sprintf('%.2f', wtau), colour = abs(wtau) > 0.5), size = 4.2, fontface = 'bold') +
  scale_fill_gradient2(low = '#b2182b', mid = '#f7f7f7', high = '#3f007d', midpoint = 0, limits = c(-1, 1), name = "Weighted Kendall's τ",
                       guide = guide_colourbar(title.position = 'top', title.hjust = 0.5)) +
  scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = 'grey15'), guide = 'none') +
  coord_equal() +
  scale_x_discrete(labels = \(x) str_replace(x, ' ', '\n')) +
  labs(x = NULL, y = NULL,
       caption = sprintf(paste0("Agreement between the tissues' rankings of\nthe %d methods (overall improvement):\n",
                                "weighted Kendall's τ (hyperbolic weights,\ndisagreements near the top count most).\n",
                                "Kendall's W = %.2f across the four tissues;\npermutation p = %.3f ",
                                "(%s permutations,\nranks shuffled within each tissue)."),
                         nrow(rmat), W, W_p, format(n_perm, big.mark = ','))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = 'grey10'), legend.position = 'bottom',
        plot.caption = element_text(hjust = 0, colour = 'grey35'), plot.caption.position = 'plot')
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('benchmark_tissue_tau', suffix, '.', ext)), pC, width = 5, height = 5, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')
# Figure 2 = strategy panel (a) + ranking consistency (b)
p2 = p1b + (pC / patchwork::plot_spacer() + patchwork::plot_layout(heights = c(1, 0.35))) +
  patchwork::plot_layout(widths = c(2.1, 1)) + patchwork::plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold', size = 14))
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('benchmark_figure2', suffix, '.', ext)), p2, width = 14.5, height = 7.5, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')

# ---- 2. appendix heatmap: scores, colour = min-max within column ----
hm_rows = scores |> filter(case == 'dissected') |>
  mutate(group = ifelse(method %in% ref_rows, 'Reference', as.character(input_of(method))))
hm = hm_rows |> select(method, group, tissue, overall, bio, batch) |>
  pivot_longer(c(overall, bio, batch), names_to = 'score', values_to = 'value')
hm = bind_rows(hm, hm |> group_by(method, group, score) |> summarise(tissue = 'mean', value = mean(value), .groups = 'drop')) |>
  group_by(score, tissue) |> mutate(scaled = (value - min(value)) / (max(value) - min(value))) |> ungroup()
grp_ord = hm |> filter(score == 'overall', tissue == 'mean') |>
  mutate(group = factor(group, levels = c(input_levels, 'Reference'))) |> arrange(group, -value)
hm = hm |> mutate(
  x = factor(ifelse(tissue == 'mean', 'Mean', str_replace(tissue_labels[tissue], ' ', '\n')),
             levels = c(str_replace(tissue_labels, ' ', '\n'), 'Mean')),
  y = factor(label_of(method), levels = rev(label_of(grp_ord$method))),
  group = factor(group, levels = c(input_levels, 'Reference')),
  score = factor(score, levels = c('overall', 'bio', 'batch'),
                 labels = c('Overall (0.6 bio + 0.4 batch)', 'Bio-conservation', 'Batch correction')))
p2 = ggplot(hm, aes(x, y)) +
  geom_tile(aes(fill = scaled), colour = 'white', linewidth = 0.6, width = 0.95, height = 0.9) +
  geom_tile(data = filter(hm, tissue == 'mean'), fill = NA, colour = 'grey20', linewidth = 0.45, width = 0.95, height = 0.9) +
  geom_text(aes(label = sprintf('%.3f', value), colour = scaled > 0.55), size = 2.5,
            fontface = ifelse(hm$tissue == 'mean', 'bold', 'plain')) +
  facet_grid(group ~ score, scales = 'free_y', space = 'free_y', switch = 'y') +
  scale_fill_gradient(low = '#f2f7fd', high = '#08306b', limits = c(0, 1), breaks = c(0, 0.5, 1),
                      labels = c('column min', '', 'column max'), name = 'Score (scaled within column)',
                      guide = guide_colourbar(title.position = 'top', title.hjust = 0.5, barwidth = unit(14, 'lines'))) +
  scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = 'grey15'), guide = 'none') +
  scale_x_discrete(position = 'top') + scale_y_discrete(position = 'right') +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid = element_blank(), panel.spacing.y = unit(0.6, 'lines'), panel.spacing.x = unit(0.8, 'lines'),
        strip.text.x = element_text(face = 'bold', size = 10.5), strip.text.y.left = element_text(face = 'bold', angle = 0, hjust = 1),
        strip.placement = 'outside', axis.text.y.right = element_text(colour = 'grey10', hjust = 0),
        legend.position = 'bottom')
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('benchmark_heatmap', suffix, '.', ext)), p2, width = 14, height = 9.5, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')

# ---- 3. orthologue removal: overall improvement vs % orthologues removed ----
if (length(missing_cases)) {
  cat('\northologue_removal', suffix, ': skipped, extra metrics missing for', paste(missing_cases, collapse = ', '), '\n')
} else {
  orth = imp |> mutate(removed = cases[case], tissue = factor(tissue, levels = names(tissue_labels), labels = tissue_labels))
  orth_mean = orth |> group_by(tissue, removed, input) |> summarise(imp = mean(imp_overall), n = n(), .groups = 'drop')
  input_line = c(Native = '#2a78d6', `scVI-latent` = '#1baf7a', UINMF = '#eb6834')
  p3 = ggplot(orth, aes(removed, imp_overall, colour = input)) +
    geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40', linewidth = 0.4) +
    geom_line(aes(group = method), linewidth = 0.3, alpha = 0.3) +
    geom_line(data = orth_mean, aes(y = imp), linewidth = 1.3) +
    geom_point(data = orth_mean, aes(y = imp), size = 2) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_x_continuous(breaks = cases, labels = \(x) paste0(x, '%')) +
    scale_colour_manual(values = input_line, name = 'Input (thick: mean over methods)') +
    labs(x = 'Shared orthologues removed (0% = dissected)', y = 'Overall improvement over unintegrated') +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = 'grey95', colour = 'grey70'),
          strip.text = element_text(face = 'bold'), legend.position = 'bottom')
  for (ext in c('pdf', 'png'))
    ggsave(file.path(out, paste0('orthologue_removal', suffix, '.', ext)), p3, width = 13, height = 4.5, dpi = 300, bg = 'white',
           device = if (ext == 'pdf') cairo_pdf else 'png')
  cat('\nOrthologue removal, methods per input:', paste(names(table(orth_mean$input)), table(orth$input[orth$case == 'dissected' & orth$tissue == 'Colon']), collapse = ', '), '\n')
}

# ---- 4. critical-difference diagram (Demsar 2006): Friedman mean ranks of overall, embedding methods, dissected ----
# Methods joined by a bar are not significantly different by the Nemenyi test (alpha = 0.05) over the four tissues.
cd_rank = dis |> filter(!method %in% graph_only) |> group_by(tissue) |> mutate(r = rank(-overall)) |>
  group_by(method, input) |> summarise(r = mean(r), .groups = 'drop') |> arrange(r)
k_m = nrow(cd_rank); n_t = n_distinct(dis$tissue)
cd = qtukey(0.95, k_m, Inf) / sqrt(2) * sqrt(k_m * (k_m + 1) / (6 * n_t))
fr = friedman.test(as.matrix(dis |> filter(!method %in% graph_only) |> select(tissue, method, overall) |>
                               pivot_wider(names_from = method, values_from = overall) |> select(-tissue)))
cat(sprintf('\nCritical difference: Friedman chi2 = %.1f, df = %d, p = %.3g; Nemenyi CD = %.1f ranks (%d methods, %d tissues)\n',
            fr$statistic, fr$parameter, fr$p.value, cd, k_m, n_t))
# maximal groups of consecutive methods spanning less than CD
grp = purrr::map(seq_len(k_m), \(i) which(cd_rank$r >= cd_rank$r[i] & cd_rank$r - cd_rank$r[i] < cd))
grp = grp[purrr::map_lgl(seq_along(grp), \(i) length(grp[[i]]) > 1 &&
                    !any(purrr::map_lgl(seq_along(grp)[-i], \(j) all(grp[[i]] %in% grp[[j]]) && length(grp[[j]]) > length(grp[[i]]))))]
grp = unique(grp)
half = ceiling(k_m / 2)
lab_df = cd_rank |> mutate(i = row_number(), side = ifelse(i <= half, 'left', 'right'),
                           y = -ifelse(side == 'left', i, k_m - i + 1) * 0.55 - 0.6 - length(grp) * 0.25,
                           xl = ifelse(side == 'left', 0.5, k_m + 0.5), label = label_of(method))
bars = purrr::imap_dfr(grp, \(g, j) tibble(x = min(cd_rank$r[g]) - 0.1, xend = max(cd_rank$r[g]) + 0.1, y = -0.25 - j * 0.25))
p4 = ggplot() +
  annotate('segment', x = 1, xend = k_m, y = 0, yend = 0, linewidth = 0.5) +
  annotate('segment', x = 1:k_m, xend = 1:k_m, y = 0, yend = 0.12, linewidth = 0.4) +
  annotate('text', x = c(1, 5, 10, 15, k_m), y = 0.32, label = c(1, 5, 10, 15, k_m), size = 3.3) +
  annotate('segment', x = 1, xend = 1 + cd, y = 0.85, yend = 0.85, linewidth = 0.9) +
  annotate('segment', x = c(1, 1 + cd), xend = c(1, 1 + cd), y = 0.78, yend = 0.92, linewidth = 0.9) +
  annotate('text', x = 1 + cd / 2, y = 1.05, label = sprintf('CD = %.1f', cd), size = 3.4) +
  geom_segment(data = lab_df, aes(x = r, xend = r, y = 0, yend = y), linewidth = 0.35, colour = 'grey30') +
  geom_segment(data = lab_df, aes(x = r, xend = xl, y = y, yend = y), linewidth = 0.35, colour = 'grey30') +
  geom_point(data = lab_df, aes(r, 0, fill = input), shape = 21, size = 2.6) +
  geom_text(data = lab_df, aes(xl, y, label = label, hjust = ifelse(side == 'left', 1.05, -0.05)), size = 3.3) +
  geom_segment(data = bars, aes(x = x, xend = xend, y = y, yend = y), linewidth = 2, colour = 'grey20', lineend = 'round') +
  scale_fill_manual(values = input_fill, name = 'Input', drop = FALSE) +
  scale_x_continuous(expand = expansion(add = c(6.5, 6.5))) +
  labs(x = 'Mean rank of the overall score over four tissues (1 = best)') +
  theme_void(base_size = 11) + theme(legend.position = 'bottom', axis.title.x = element_text(margin = margin(t = 6)))
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('critical_difference', suffix, '.', ext)), p4, width = 12, height = 6, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')
cat('CD groups (methods not separable):', length(grp), '\n')
