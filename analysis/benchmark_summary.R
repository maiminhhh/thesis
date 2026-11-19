# Benchmark summary across all tissues and cases:  Rscript analysis/benchmark_summary.R
# 4 tissues x 6 cases (full, dissected, genes2-8) = 24 datasets, 23 integration methods.
# Overall = 0.6 bio + 0.4 batch; bio uses the four metrics every tissue has (NMI, ARI,
# ASW_label, cLISI). Methods are ranked within each dataset (1 = best), so tissues with
# different score levels count equally.
# Output: figures/benchmark_summary/*.png and scores_all.csv
source('library.R')
library(ggrepel)
library(patchwork)

out = 'figures/benchmark_summary'
dir.create(out, showWarnings = FALSE, recursive = TRUE)

tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon',
                  breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
case_levels = c('full', 'dissected', 'genes2', 'genes4', 'genes6', 'genes8')
score_dir = function(tissue, case) {
  if (tissue == 'breast_cancer')
    return(sprintf('breast_cancer/data/%s/', if (case %in% c('full', 'dissected')) paste0(case, '_mouse') else case))
  if (case %in% c('full', 'dissected')) sprintf('%s/data/05_evaluated_score/%s_mouse/', tissue, case)
  else sprintf('%s/data/%s/', tissue, case)
}
# shared one-to-one orthologues per dataset (full uses the same genes as dissected)
n_genes = tribble(~tissue, ~dissected, ~genes2, ~genes4, ~genes6, ~genes8,
                  'bone_marrow',   123, 98, 74, 49, 25,
                  'colon',          46, 37, 28, 18,  9,
                  'breast_cancer', 122, 98, 73, 49, 24,
                  'fetal_liver',   123, 98, 74, 49, 25) |>
  mutate(full = dissected) |>
  pivot_longer(-tissue, names_to = 'case', values_to = 'n_shared')

batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI')
read_metrics = function(f) {
  df = suppressMessages(read_csv(f, show_col_types = FALSE))
  tibble(metric = df[[1]], value = as.numeric(df[[2]]), method = names(df)[2])
}

scores = expand_grid(tissue = names(tissue_labels), case = case_levels) |>
  pmap_dfr(function(tissue, case) {
    list.files(score_dir(tissue, case), pattern = '^05\\..*\\.csv$', full.names = TRUE) |>
      map_dfr(read_metrics) |>
      mutate(group = case_when(metric %in% batch_metrics ~ 'batch',
                               metric %in% bio_metrics   ~ 'bio')) |>
      filter(!is.na(group), !is.na(value)) |>
      group_by(method, group) |> summarise(score = mean(value), .groups = 'drop') |>
      pivot_wider(names_from = group, values_from = score) |>
      mutate(tissue = tissue, case = case)
  }) |>
  mutate(overall = 0.6 * bio + 0.4 * batch,
         input = case_when(str_detect(method, 'UINMF')       ~ 'UINMF',
                           str_detect(method, 'scVI-latent') ~ 'scVI latent',
                           TRUE                              ~ 'PCA / own'),
         input = factor(input, levels = c('PCA / own', 'UINMF', 'scVI latent')),
         base  = str_remove(method, ' (UINMF|scVI-latent)$')) |>
  left_join(n_genes, by = c('tissue', 'case'))

methods = scores |> filter(!str_starts(method, 'Unintegrated')) |>
  group_by(tissue, case) |>
  mutate(rank = rank(-overall), rank_bio = rank(-bio), rank_batch = rank(-batch)) |> ungroup()
stopifnot(nrow(methods) == 24 * 23)
write_csv(methods, file.path(out, 'scores_all.csv'))

summary = methods |> group_by(method, input) |>
  summarise(mean_rank = mean(rank), median_rank = median(rank),
            se_bio = sd(bio) / sqrt(n()), se_batch = sd(batch) / sqrt(n()),
            overall = mean(overall), bio = mean(bio), batch = mean(batch), .groups = 'drop') |>
  arrange(mean_rank)
method_order = rev(summary$method)

# validated categorical palette (slots 1-3 all-pairs safe; slot 4 needs the shape encoding)
input_cols    = c('PCA / own' = '#2a78d6', 'UINMF' = '#eb6834', 'scVI latent' = '#1baf7a')
input_shapes  = c('PCA / own' = 16, 'UINMF' = 17, 'scVI latent' = 15)
tissue_cols   = setNames(c('#2a78d6', '#eb6834', '#1baf7a', '#eda100'), tissue_labels)
tissue_shapes = setNames(c(16, 17, 15, 18), tissue_labels)
blue_ramp     = c('#104281', '#256abf', '#5598e7', '#9ec5f4', '#e8f1fc')   # dark = good rank

theme_bench = function(base_size = 11) theme_bw(base_size = base_size) +
  theme(panel.grid.minor = element_blank(), plot.title.position = 'plot',
        strip.background = element_rect(fill = 'grey95', colour = 'grey70'),
        strip.text = element_text(face = 'bold'), legend.position = 'top')
save = function(p, file, w, h) ggsave(file.path(out, file), p, width = w, height = h, dpi = 300, bg = 'white')

# ---- 1. summary heatmap: rank per tissue + mean over datasets, with score bars ----
# m: the method rows to summarise (one or more datasets per tissue)
summary_heatmap = function(m, all_label, title, subtitle, file) {
  summ = m |> group_by(method, input) |>
    summarise(mean_rank = mean(rank), overall = mean(overall), bio = mean(bio), batch = mean(batch),
              .groups = 'drop') |> arrange(mean_rank)
  ord = rev(summ$method)
  tr = m |> group_by(method, tissue) |> summarise(r = mean(rank), .groups = 'drop') |>
    mutate(column = tissue_labels[tissue]) |>
    bind_rows(summ |> transmute(method, r = mean_rank, column = all_label)) |>
    mutate(column = factor(column, levels = c(tissue_labels, all_label)),
           method = factor(method, levels = ord),
           lab = if_else(r == round(r), sprintf('%d', as.integer(r)), sprintf('%.1f', r)))
  h1 = ggplot(tr, aes(column, method, fill = r)) +
    geom_tile(colour = 'white', linewidth = 0.8) +
    geom_text(aes(label = lab, colour = r < 9), size = 3.2) +
    scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = '#0b0b0b'), guide = 'none') +
    scale_fill_gradientn(colours = blue_ramp, limits = c(1, 23), name = 'Rank\n(1 = best)') +
    scale_x_discrete(position = 'top') +
    labs(x = NULL, y = NULL) + theme_bench() +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(), legend.position = 'bottom',
          axis.text.x.top = element_text(face = 'bold'))
  h0 = ggplot(summ |> mutate(method = factor(method, levels = ord)),
              aes(x = 1, y = method, colour = input, shape = input)) +
    geom_point(size = 3) +
    scale_colour_manual(values = input_cols) + scale_shape_manual(values = input_shapes) +
    theme_void() + theme(legend.position = 'none')
  bars = summ |> select(method, Batch = batch, Bio = bio, Overall = overall) |>
    pivot_longer(-method, names_to = 'metric') |>
    mutate(method = factor(method, levels = ord), metric = factor(metric, levels = c('Batch', 'Bio', 'Overall')))
  h2 = ggplot(bars, aes(value, method)) +
    geom_col(fill = '#5598e7', width = 0.65) +
    geom_text(aes(label = sprintf('%.2f', value)), hjust = -0.15, size = 2.9, colour = '#52514e') +
    facet_grid(. ~ metric) + scale_x_continuous(limits = c(0, 0.85), breaks = c(0, 0.4, 0.8)) +
    labs(x = if (all_label == 'All 24') 'Mean score over 24 datasets' else 'Mean score over the four tissues', y = NULL) + theme_bench() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
  p = (h0 + h1 + h2 + plot_layout(widths = c(0.04, 1, 1.1))) +
    plot_annotation(title = title, subtitle = subtitle, theme = theme(plot.title = element_text(face = 'bold')))
  save(p, file, 15, 9)
  summ
}
dot_note = 'Dot = input type: blue PCA / own, orange UINMF, green scVI latent.'
summ_dissected = summary_heatmap(
  methods |> filter(case == 'dissected'), 'Mean of 4',
  'Integration benchmark summary: 23 methods on the four dissected tissues',
  paste('Rank within each tissue by overall score (0.6 bio + 0.4 batch); last column: mean rank over the four tissues.', dot_note),
  '1_summary_heatmap.png')
summary_heatmap(methods, 'All 24',
  'Integration benchmark summary: 23 methods, 4 tissues x 6 datasets',
  paste('Rank within each dataset by overall score (0.6 bio + 0.4 batch), averaged per tissue and over all 24 datasets.', dot_note),
  '1b_summary_heatmap_all24.png') |> invisible()

# ---- 2. rank distribution over the 24 datasets ----
p2 = ggplot(methods |> mutate(method = factor(method, levels = method_order)), aes(rank, method)) +
  geom_boxplot(aes(fill = input), alpha = 0.25, outlier.shape = NA, width = 0.65, colour = 'grey40') +
  geom_point(aes(colour = factor(tissue_labels[tissue], levels = tissue_labels), shape = factor(tissue_labels[tissue], levels = tissue_labels)), size = 1.6, alpha = 0.8,
             position = position_jitter(height = 0.18, width = 0, seed = 1)) +
  scale_fill_manual(values = input_cols, name = 'Input') +
  scale_colour_manual(values = tissue_cols, name = 'Tissue') + scale_shape_manual(values = tissue_shapes, name = 'Tissue') +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 23)) +
  labs(x = 'Rank within dataset (1 = best)', y = NULL,
       title = 'How consistent is each method? Rank over 24 datasets',
       subtitle = 'Narrow box near 1: reliably good. Wide box: depends on the dataset.') +
  theme_bench() + theme(panel.grid.major.y = element_blank())
save(p2, '2_rank_distribution.png', 11, 9)

# ---- 3. robustness to fewer shared genes (same six base methods on each input) ----
six = c('FastMNN', 'CCA', 'Harmony', 'Scanorama', 'RPCA', 'BBKNN')
gene_d = methods |> filter(base %in% six, case != 'full') |>
  group_by(tissue, case, n_shared, input) |>
  summarise(overall = mean(overall), lo = min(overall), .groups = 'drop')
gene_ref = methods |> filter(method %in% c('LIGER UINMF', 'scVI', 'SAMap'), case != 'full')
p3 = ggplot(gene_d, aes(n_shared, overall, colour = input, shape = input)) +
  geom_line(linewidth = 1) + geom_point(size = 2.6) +
  facet_wrap(~ factor(tissue_labels[tissue], levels = tissue_labels), nrow = 1, scales = 'free_x') +
  scale_colour_manual(values = input_cols, name = 'Input (mean of FastMNN, CCA, Harmony, Scanorama, RPCA, BBKNN)') +
  scale_shape_manual(values = input_shapes, name = 'Input (mean of FastMNN, CCA, Harmony, Scanorama, RPCA, BBKNN)') +
  scale_x_reverse() +
  labs(x = 'Shared one-to-one orthologues (dissected, then 20-80% removed)', y = 'Overall score',
       title = 'Does performance hold up with fewer shared genes?',
       subtitle = 'Same six integration methods run on each input; x axis runs from most to fewest shared genes.') +
  theme_bench()
save(p3, '3_gene_removal.png', 13, 4.8)

# ---- 4. batch vs bio trade-off ----
p4 = ggplot(summary, aes(batch, bio, colour = input, shape = input)) +
  geom_errorbar(aes(ymin = bio - se_bio, ymax = bio + se_bio), width = 0, alpha = 0.5) +
  geom_errorbar(aes(xmin = batch - se_batch, xmax = batch + se_batch), width = 0, alpha = 0.5, orientation = 'y') +
  geom_point(size = 3) +
  geom_text_repel(aes(label = method), size = 3, colour = '#0b0b0b', seed = 1, max.overlaps = 30, box.padding = 0.4) +
  scale_colour_manual(values = input_cols, name = 'Input') + scale_shape_manual(values = input_shapes, name = 'Input') +
  labs(x = 'Batch removal (mean over 24 datasets)', y = 'Bio conservation (mean over 24 datasets)',
       title = 'Trade-off between removing species effects and keeping cell types apart',
       subtitle = 'Error bars: +/- 1 standard error over datasets. Top right is best.') +
  theme_bench()
save(p4, '4_batch_vs_bio.png', 10, 8)

# ---- 5. does the input help? paired difference vs the same method on PCA ----
paired = methods |> filter(base %in% six) |>
  select(tissue, case, base, input, overall) |>
  pivot_wider(names_from = input, values_from = overall) |>
  pivot_longer(c(UINMF, `scVI latent`), names_to = 'input', values_to = 'alt') |>
  mutate(delta = alt - `PCA / own`, input = factor(input, levels = c('UINMF', 'scVI latent')),
         tissue = factor(tissue_labels[tissue], levels = tissue_labels))
p5 = ggplot(paired, aes(delta, fct_reorder(base, delta))) +
  geom_vline(xintercept = 0, colour = 'grey55') +
  geom_boxplot(fill = 'grey92', colour = 'grey45', outlier.shape = NA, width = 0.6) +
  geom_point(aes(colour = tissue, shape = tissue), size = 1.8, alpha = 0.85,
             position = position_jitter(height = 0.18, width = 0, seed = 1)) +
  facet_grid(. ~ input) +
  scale_colour_manual(values = tissue_cols, name = 'Tissue') + scale_shape_manual(values = tissue_shapes, name = 'Tissue') +
  labs(x = 'Overall score: method on this input minus the same method on PCA', y = NULL,
       title = 'Does changing the input help each method?',
       subtitle = 'One point per dataset (24). Right of the line: the alternative input beats PCA.') +
  theme_bench() + theme(panel.grid.major.y = element_blank())
save(p5, '5_input_effect.png', 11, 5.5)

# ---- 6. UINMF advantage vs number of shared genes ----
gain = paired |> group_by(tissue, case, input) |> summarise(delta = mean(delta), .groups = 'drop') |>
  left_join(n_genes |> mutate(tissue = factor(tissue_labels[tissue], levels = tissue_labels)), by = c('tissue', 'case'))
p6 = ggplot(gain, aes(n_shared, delta)) +
  geom_hline(yintercept = 0, colour = 'grey55') +
  geom_smooth(method = 'lm', formula = y ~ x, colour = 'grey40', fill = 'grey85', linewidth = 0.7) +
  geom_point(aes(colour = tissue, shape = tissue), size = 2.8) +
  facet_grid(. ~ input) +
  scale_colour_manual(values = tissue_cols, name = 'Tissue') + scale_shape_manual(values = tissue_shapes, name = 'Tissue') +
  labs(x = 'Shared one-to-one orthologues', y = 'Mean gain over PCA (six methods)',
       title = 'Is the alternative input more useful when fewer genes are shared?',
       subtitle = 'One point per dataset (24); line: linear fit with 95% band.') +
  theme_bench()
save(p6, '6_gain_vs_shared_genes.png', 11, 5)

print(summary |> mutate(across(where(is.numeric), \(x) round(x, 3))), n = Inf)
print(gain |> group_by(input) |> summarise(cor = cor(n_shared, delta, method = 'spearman')))

# ---- 7. absolute overall score per method over the 24 datasets ----
# the three unintegrated inputs are shown as reference rows (grey)
abs_d = scores |>
  mutate(kind = if_else(str_starts(method, 'Unintegrated'), 'Unintegrated', as.character(input)),
         kind = factor(kind, levels = c('PCA / own', 'UINMF', 'scVI latent', 'Unintegrated')),
         tissue = factor(tissue_labels[tissue], levels = tissue_labels)) |>
  group_by(method) |> mutate(med = median(overall)) |> ungroup() |>
  mutate(method = fct_reorder(method, med))
p7 = ggplot(abs_d, aes(overall, method)) +
  geom_vline(xintercept = median(abs_d$overall[abs_d$method == 'Unintegrated']),
             linetype = 'dashed', colour = 'grey50') +
  geom_boxplot(aes(fill = kind), alpha = 0.3, outlier.shape = NA, width = 0.65, colour = 'grey40') +
  geom_point(aes(colour = tissue, shape = tissue), size = 1.7, alpha = 0.85,
             position = position_jitter(height = 0.18, width = 0, seed = 1)) +
  scale_fill_manual(values = c(input_cols, Unintegrated = '#9a9a96'), name = 'Input') +
  scale_colour_manual(values = tissue_cols, name = 'Tissue') + scale_shape_manual(values = tissue_shapes, name = 'Tissue') +
  labs(x = 'Overall score (0.6 bio + 0.4 batch)', y = NULL,
       title = 'Overall integration score by method, over 24 datasets',
       subtitle = 'One point per dataset (4 tissues x 6 datasets). Dashed line: median of unintegrated PCA.') +
  theme_bench() + theme(panel.grid.major.y = element_blank(), legend.box = 'vertical')
save(p7, '7_overall_score_boxplot.png', 11, 10)

# ---- 8. UINMF vs scVI-based models as shared genes are removed (dissected chain) ----
level_pct = c(dissected = 100, genes2 = 80, genes4 = 60, genes6 = 40, genes8 = 20)
facet_lab = n_genes |> filter(case %in% c('dissected', 'genes8')) |>
  pivot_wider(names_from = case, values_from = n_shared) |>
  transmute(tissue, lab = sprintf('%s (%d to %d genes)', tissue_labels[tissue], dissected, genes8))
facet_lab = setNames(facet_lab$lab, facet_lab$tissue)
chain = methods |> filter(case != 'full') |> mutate(pct = level_pct[case],
  tissue_f = factor(facet_lab[tissue], levels = facet_lab[names(tissue_labels)]))

# A: same six integration methods, UINMF input minus scVI-latent input
pair6 = chain |> filter(base %in% six, input != 'PCA / own') |>
  select(tissue_f, case, pct, base, input, overall, bio, batch) |>
  pivot_wider(names_from = input, values_from = c(overall, bio, batch)) |>
  mutate(d_overall = overall_UINMF - `overall_scVI latent`,
         d_bio = bio_UINMF - `bio_scVI latent`, d_batch = batch_UINMF - `batch_scVI latent`)
pA = ggplot(pair6, aes(pct, d_overall)) +
  geom_hline(yintercept = 0, colour = 'grey55') +
  geom_point(colour = '#52514e', alpha = 0.55, size = 1.8, position = position_jitter(width = 2, height = 0, seed = 1)) +
  stat_summary(fun = mean, geom = 'line', colour = '#eb6834', linewidth = 1.1) +
  stat_summary(fun = mean, geom = 'point', colour = '#eb6834', size = 2.8) +
  facet_wrap(~ tissue_f, nrow = 1) +
  scale_x_reverse(breaks = c(100, 80, 60, 40, 20)) +
  labs(x = NULL, y = 'UINMF minus scVI latent\n(overall score)',
       title = 'A. Same six methods (FastMNN, CCA, Harmony, Scanorama, RPCA, BBKNN) on each input',
       subtitle = 'Grey points: one method; orange line: mean. Above 0: UINMF input better.') +
  theme_bench()

# B: the standalone models, absolute overall score
flag_cols   = c('LIGER UINMF' = '#eb6834', 'scVI' = '#1baf7a', 'scANVI' = '#2a78d6')
flag_shapes = c('LIGER UINMF' = 17, 'scVI' = 15, 'scANVI' = 16)
pB = ggplot(chain |> filter(method %in% names(flag_cols)), aes(pct, overall, colour = method, shape = method)) +
  geom_line(linewidth = 1) + geom_point(size = 2.6) +
  facet_wrap(~ tissue_f, nrow = 1) +
  scale_colour_manual(values = flag_cols, name = NULL) + scale_shape_manual(values = flag_shapes, name = NULL) +
  scale_x_reverse(breaks = c(100, 80, 60, 40, 20)) +
  labs(x = 'Shared one-to-one orthologues kept (%; dissected = 100)', y = 'Overall score',
       title = 'B. Standalone models: LIGER UINMF vs scVI and scANVI') +
  theme_bench()
save(pA / pB, '8_uinmf_vs_scvi_gene_removal.png', 14, 9)

# C: which component drives the difference (six-method mean)
compC = pair6 |> group_by(tissue_f, pct) |>
  summarise(Batch = mean(d_batch), Bio = mean(d_bio), .groups = 'drop') |>
  pivot_longer(c(Batch, Bio), names_to = 'component', values_to = 'delta')
pC = ggplot(compC, aes(pct, delta, colour = component, shape = component)) +
  geom_hline(yintercept = 0, colour = 'grey55') +
  geom_line(linewidth = 1) + geom_point(size = 2.6) +
  facet_wrap(~ tissue_f, nrow = 1) +
  scale_colour_manual(values = c(Batch = '#2a78d6', Bio = '#eb6834'), name = NULL) +
  scale_shape_manual(values = c(Batch = 16, Bio = 17), name = NULL) +
  scale_x_reverse(breaks = c(100, 80, 60, 40, 20)) +
  labs(x = 'Shared one-to-one orthologues kept (%; dissected = 100)', y = 'UINMF minus scVI latent',
       title = 'Which part of the score differs? Six-method mean, batch and bio components',
       subtitle = 'Above 0: UINMF input better on that component.') +
  theme_bench()
save(pC, '8b_uinmf_vs_scvi_components.png', 14, 4.8)

win = pair6 |> group_by(tissue_f, pct) |>
  summarise(mean_diff = round(mean(d_overall), 3), uinmf_wins = sum(d_overall > 0), of = n(), .groups = 'drop')
write_csv(win, file.path(out, 'uinmf_vs_scvi_latent.csv'))
print(win, n = Inf)
print(pair6 |> group_by(pct) |> summarise(mean_diff = round(mean(d_overall), 3), uinmf_wins = sum(d_overall > 0), of = n()))
print(chain |> filter(method %in% names(flag_cols)) |> select(tissue, pct, method, overall) |>
  pivot_wider(names_from = method, values_from = overall) |> mutate(across(where(is.double), \(x) round(x, 3))), n = Inf)

# ---- 9. BENGAL-style figure: (a) integrated score per strategy, (b) mixing vs bio per tissue ----
# colour = integration method (base, colours as in analysis/visualisations.R), shape = input
base_cols9 = c(Unintegrated = '#000000', FastMNN = '#E69F00', SAMap = '#0072B2', CCA = '#009E73',
               Scanorama = '#D55E00', RPCA = '#332288', Harmony = '#AA4499', BBKNN = '#999933',
               scVI = '#882255', scANVI = '#CC6677', LIGER = '#8C510A')
input_shapes9 = c('PCA / own' = 15, 'UINMF' = 17, 'scVI latent' = 3)
strat = scores |>
  mutate(base = case_when(str_starts(method, 'Unintegrated') ~ 'Unintegrated',
                          str_starts(method, 'LIGER')        ~ 'LIGER',
                          TRUE ~ base),
         input = case_when(method == 'Unintegrated UINMF' ~ 'UINMF',
                           method == 'Unintegrated scVI-latent' ~ 'scVI latent',
                           TRUE ~ as.character(input)),
         input = factor(input, levels = names(input_shapes9)),
         base = factor(base, levels = names(base_cols9)),
         label = case_when(method == 'Unintegrated' ~ 'Unintegrated PCA',
                           TRUE ~ str_replace(method, 'scVI-latent', 'scVI latent')))

# gain = overall score minus that of the raw unintegrated data (PCA of shared genes) in the
# same dataset, for every strategy; the unintegrated UINMF and scVI-latent inputs are kept,
# so their bars show what the input alone adds
raw = strat |> filter(method == 'Unintegrated') |> select(tissue, case, overall0 = overall)
a_d = strat |> filter(method != 'Unintegrated') |>
  left_join(raw, by = c('tissue', 'case')) |>
  mutate(gain = overall - overall0) |>
  group_by(label) |> mutate(med = median(gain)) |> ungroup() |>
  mutate(label = fct_reorder(label, -med))
stopifnot(!anyNA(a_d$gain))
# panel a is coloured by what each strategy is built on
built_cols = c(PCA = '#2a78d6', UINMF = '#eb6834', scVI = '#1baf7a', 'Other (LIGER iNMF, SAMap)' = '#eda100')
a_d = a_d |> mutate(built = case_when(
  str_detect(method, 'UINMF')                                     ~ 'UINMF',
  str_detect(method, 'scVI') | method == 'scANVI'                  ~ 'scVI',
  method %in% c('LIGER', 'SAMap')                                  ~ 'Other (LIGER iNMF, SAMap)',
  TRUE                                                             ~ 'PCA'),
  built = factor(built, levels = names(built_cols)))
p9a = ggplot(a_d, aes(label, gain, fill = built)) +
  geom_hline(yintercept = 0, colour = 'grey50', linetype = 'dashed') +
  geom_boxplot(colour = 'grey20', linewidth = 0.35, width = 0.75, outlier.size = 1) +
  scale_fill_manual(values = built_cols, name = 'Built on') +
  labs(x = 'Integration strategy', y = 'Gain over unintegrated raw data\n(overall score, 0.6 bio + 0.4 batch)', tag = 'a') +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9.5),
        legend.position = 'top', plot.tag = element_text(face = 'bold', size = 16))

# Pareto front per panel: strategies no other (integrated) strategy beats on both batch and bio
pareto_front = function(x) {
  x = x[order(-x$batch, -x$bio), ]
  x[x$bio > c(-Inf, head(cummax(x$bio), -1)), ]
}
scatter_panel = function(d, facet_var, tag, group_vars = facet_var, pareto = FALSE) {
  front = d |> filter(base != 'Unintegrated') |> group_by(across(all_of(group_vars))) |>
    group_modify(~ pareto_front(.x)) |> ungroup() |> arrange(batch)
  med = d |> group_by(across(all_of(group_vars))) |>
    summarise(un_b = median(batch[base == 'Unintegrated']), un_o = median(bio[base == 'Unintegrated']),
              all_b = median(batch), all_o = median(bio), .groups = 'drop')
  ggplot(d, aes(batch, bio, colour = base, shape = input)) +
    geom_vline(data = med, aes(xintercept = un_b), linetype = 'dashed', colour = '#d8a1b5', linewidth = 0.6) +
    geom_hline(data = med, aes(yintercept = un_o), linetype = 'dashed', colour = '#d8a1b5', linewidth = 0.6) +
    geom_vline(data = med, aes(xintercept = all_b), linetype = 'dashed', colour = 'grey60', linewidth = 0.6) +
    geom_hline(data = med, aes(yintercept = all_o), linetype = 'dashed', colour = 'grey60', linewidth = 0.6) +
    { if (pareto) geom_step(data = front, aes(batch, bio), direction = 'hv', inherit.aes = FALSE,
                            colour = 'grey35', linewidth = 0.5) } +
    geom_point(size = 2.8, stroke = 1.1, alpha = 0.9) +
    { if (pareto) geom_text_repel(data = front, aes(batch, bio, label = label), inherit.aes = FALSE,
                                  size = 3, colour = '#0b0b0b', seed = 1, box.padding = 0.5,
                                  min.segment.length = 0, segment.colour = 'grey50', max.overlaps = Inf) } +
    facet_wrap(as.formula(paste('~', facet_var)), nrow = 1, scales = 'free') +
    scale_colour_manual(values = base_cols9, name = 'Integration method', drop = FALSE) +
    scale_shape_manual(values = input_shapes9, name = 'Input') +
    labs(x = 'Species mixing score (batch)', y = 'Biology conservation score', tag = tag,
         caption = paste0('Dashed pink: median of the unintegrated inputs; dashed grey: median of all strategies.',
                     if (pareto) ' Solid grey: Pareto front (labelled), strategies not beaten on both scores by any other.' else '')) +
    theme_classic(base_size = 12) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 13),
          legend.position = 'right', plot.tag = element_text(face = 'bold', size = 16), panel.spacing.x = unit(1.2, 'lines'),
          plot.caption = element_text(hjust = 0, colour = 'grey35'))
}
b_d = strat |> filter(case == 'dissected') |>
  mutate(tissue = factor(tissue_labels[tissue], levels = tissue_labels))
p9b = scatter_panel(b_d, 'tissue', 'b', pareto = TRUE) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(shape = 16, size = 3.5)),
         shape = guide_legend(override.aes = list(size = 3.5)))
save(p9a / p9b + plot_layout(heights = c(1.15, 1)), '9_bengal_style.png', 16, 13)

# supplementary: (b) for every tissue x dataset
s_d = strat |> mutate(tissue = factor(tissue_labels[tissue], levels = tissue_labels),
                      case = factor(case, levels = case_levels))
p9s = scatter_panel(s_d, 'case', NULL, group_vars = c('tissue', 'case')) + facet_grid(tissue ~ case) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.3, 0.7)) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(shape = 16, size = 3.5)))
save(p9s, '9s_mixing_vs_bio_all_datasets.png', 18, 12)

# ---- 10. summary table (bars + per-metric circles), bio vs batch scatter, bio/batch boxplots ----
# dissected datasets only; per method, each metric is averaged over the four tissues
diss_raw = expand_grid(tissue = names(tissue_labels)) |>
  pmap_dfr(function(tissue) list.files(score_dir(tissue, 'dissected'), pattern = '^05\\..*\\.csv$', full.names = TRUE) |>
             map_dfr(read_metrics) |> mutate(tissue = tissue)) |>
  filter(!str_starts(method, 'Unintegrated'))
metric_cols = tribble(~key, ~group, ~label,
  'PCR_batch', 'batch', 'PCR', 'ASW_label/batch', 'batch', 'Batch ASW', 'iLISI', 'batch', 'iLISI',
  'graph_conn', 'batch', 'Graph conn.',
  'NMI_cluster/label', 'bio', 'NMI', 'ARI_cluster/label', 'bio', 'ARI', 'ASW_label', 'bio', 'Cell-type ASW',
  'cLISI', 'bio', 'cLISI')
per_metric = diss_raw |> filter(metric %in% metric_cols$key) |>
  group_by(method, metric) |> summarise(value = mean(value, na.rm = TRUE), .groups = 'drop') |>
  mutate(value = if_else(is.nan(value), NA_real_, value))
per_sum = methods |> filter(case == 'dissected') |> group_by(method, input) |>
  summarise(Overall = mean(overall), Batch = mean(batch), Bio = mean(bio), .groups = 'drop') |>
  arrange(desc(Overall)) |> mutate(y = row_number())
n_m = nrow(per_sum)

ramp = list(overall = c('#ecebe8', '#3d3c39'), batch = c('#e3eefc', '#104281'), bio = c('#fde6dc', '#9c3510'))
rank_fill = function(v, group) {                  # darker = better rank among the methods
  r = rank(-v, na.last = 'keep', ties.method = 'average')
  pal = colorRampPalette(ramp[[group]])(n_m)
  ifelse(is.na(r), NA, pal[n_m + 1 - round(r)])
}
bar_w = 1.5
layout_x = tribble(~key, ~group, ~label, ~x, ~kind,
  'Overall', 'overall', 'Overall', 0, 'bar',
  'Batch', 'batch', 'Batch', 2.1, 'bar',
  'PCR_batch', 'batch', 'PCR', 4.2, 'dot', 'ASW_label/batch', 'batch', 'Batch ASW', 5.0, 'dot',
  'iLISI', 'batch', 'iLISI', 5.8, 'dot', 'graph_conn', 'batch', 'Graph conn.', 6.6, 'dot',
  'Bio', 'bio', 'Bio', 7.6, 'bar',
  'NMI_cluster/label', 'bio', 'NMI', 9.7, 'dot', 'ARI_cluster/label', 'bio', 'ARI', 10.5, 'dot',
  'ASW_label', 'bio', 'Cell-type ASW', 11.3, 'dot', 'cLISI', 'bio', 'cLISI', 12.1, 'dot')

bars10 = per_sum |> pivot_longer(c(Overall, Batch, Bio), names_to = 'key', values_to = 'value') |>
  left_join(layout_x, by = 'key') |> group_by(key) |> mutate(fill = rank_fill(value, group[1])) |> ungroup()
dots10 = per_metric |> rename(key = metric) |> inner_join(per_sum |> select(method, y), by = 'method') |>
  left_join(layout_x, by = 'key') |> group_by(key) |> mutate(fill = rank_fill(value, group[1])) |> ungroup()
missing10 = dots10 |> filter(is.na(value))
heads = layout_x |> mutate(xc = if_else(kind == 'bar', x + bar_w / 2, x))
groups10 = tibble(label = c('Overall', 'Batch correction', 'Bio conservation'),
                  x0 = c(-0.1, 2.0, 7.5), x1 = c(1.7, 7.0, 12.5))
input_mark = per_sum |> mutate(col = input_cols[as.character(input)], shp = input_shapes[as.character(input)])
# rank legend: one small ramp per group, at the right
leg = expand_grid(group = c('overall', 'batch', 'bio'), step = 1:5) |>
  mutate(gx = 13.4 + match(group, c('overall', 'batch', 'bio')) * 0.38,
         fill = map2_chr(group, step, \(g, s) colorRampPalette(ramp[[g]])(5)[6 - s]), gy = 1.6 + step * 0.75)

p10a = ggplot() +
  geom_segment(data = groups10, aes(x = x0, xend = x1, y = -0.55, yend = -0.55), linewidth = 0.5) +
  geom_text(data = groups10, aes(x = (x0 + x1) / 2, y = -1.35, label = label), fontface = 'bold', size = 4.2) +
  geom_text(data = heads, aes(x = xc, y = 0, label = label), size = 3.1, angle = 0, colour = '#52514e') +
  geom_rect(data = bars10, aes(xmin = x, xmax = x + bar_w * value, ymin = y - 0.33, ymax = y + 0.33, fill = fill),
            colour = 'grey35', linewidth = 0.25) +
  geom_text(data = bars10, aes(x = x + bar_w * value + 0.06, y = y, label = sprintf('%.2f', value)),
            hjust = 0, size = 2.6, colour = '#52514e') +
  geom_point(data = dots10 |> filter(!is.na(value)), aes(x = x, y = y, size = value, fill = fill),
             shape = 21, colour = 'grey30', stroke = 0.3) +
  geom_text(data = missing10, aes(x = x, y = y), label = '–', colour = 'grey60', size = 3) +
  geom_point(data = input_mark, aes(x = -0.45, y = y, colour = col, shape = shp), size = 2.4) +
  geom_tile(data = leg, aes(gx, gy, fill = fill), width = 0.36, height = 0.75, colour = 'white') +
  annotate('text', x = 13.4 + 0.76, y = 0.9, label = 'Rank', fontface = 'bold', size = 3.3) +
  annotate('text', x = 14.85, y = 2.35, label = '1 (best)', hjust = 0, size = 2.8) +
  annotate('text', x = 14.85, y = 5.35, label = 'last', hjust = 0, size = 2.8) +
  scale_fill_identity() + scale_colour_identity() + scale_shape_identity() +
  scale_size_area(max_size = 7, limits = c(0, 1), breaks = c(0.2, 0.4, 0.6, 0.8, 1), name = 'Score') +
  scale_y_reverse(breaks = per_sum$y, labels = per_sum$method, expand = expansion(add = c(0.7, 0.3))) +
  scale_x_continuous(limits = c(-0.6, 15.6), expand = c(0, 0)) +
  guides(size = guide_legend(override.aes = list(fill = 'grey80'), title.position = 'top')) +
  coord_cartesian(clip = 'off') +
  labs(x = NULL, y = NULL, tag = 'a') +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = '#0b0b0b'),
        legend.position = c(0.965, 0.4), plot.tag = element_text(face = 'bold', size = 16),
        plot.margin = margin(30, 10, 5, 5))

# b: bio vs batch, mean +/- SE over the four tissues
sc = methods |> filter(case == 'dissected') |>
  mutate(bcol = case_when(str_starts(method, 'LIGER') ~ 'LIGER', TRUE ~ base)) |>
  group_by(method, input, bcol) |>
  summarise(mb = mean(batch), sb = sd(batch) / sqrt(n()), mo = mean(bio), so = sd(bio) / sqrt(n()), .groups = 'drop')
raw10 = scores |> filter(case == 'dissected', method == 'Unintegrated')
p10b = ggplot(sc, aes(mb, mo, colour = bcol, shape = input)) +
  geom_vline(xintercept = median(raw10$batch), linetype = 'dashed', colour = '#2a78d6') +
  geom_hline(yintercept = median(raw10$bio), linetype = 'dashed', colour = '#eb6834') +
  geom_errorbar(aes(ymin = mo - so, ymax = mo + so), width = 0, alpha = 0.6) +
  geom_errorbar(aes(xmin = mb - sb, xmax = mb + sb), width = 0, alpha = 0.6, orientation = 'y') +
  geom_point(size = 2.8, stroke = 1.1) +
  geom_text_repel(aes(label = method), size = 2.6, colour = '#0b0b0b', seed = 1, max.overlaps = Inf,
                  box.padding = 0.35, min.segment.length = 0, segment.colour = 'grey60') +
  scale_colour_manual(values = base_cols9, name = 'Integration method') +
  scale_shape_manual(values = input_shapes9, name = 'Input') +
  labs(x = 'Batch correction (mean over 4 tissues)', y = 'Bio conservation (mean over 4 tissues)', tag = 'b',
       caption = 'Bars: +/- 1 SE over tissues. Dashed: unintegrated data (median over tissues).') +
  theme_classic(base_size = 11) +
  theme(legend.position = 'bottom', legend.box = 'vertical', plot.tag = element_text(face = 'bold', size = 16),
        plot.caption = element_text(hjust = 0, colour = 'grey35')) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(shape = 16)))

# c: bio and batch per method over the four tissues
bx = methods |> filter(case == 'dissected') |> select(method, tissue, Batch = batch, Bio = bio) |>
  pivot_longer(c(Batch, Bio), names_to = 'component', values_to = 'score') |>
  mutate(method = factor(method, levels = per_sum$method))
p10c = ggplot(bx, aes(method, score, fill = component, colour = component)) +
  geom_hline(yintercept = median(raw10$batch), linetype = 'dashed', colour = '#2a78d6') +
  geom_hline(yintercept = median(raw10$bio), linetype = 'dashed', colour = '#eb6834') +
  geom_boxplot(alpha = 0.2, width = 0.7, outlier.size = 0.8, linewidth = 0.4,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c(Batch = '#2a78d6', Bio = '#eb6834'), name = NULL) +
  scale_colour_manual(values = c(Batch = '#2a78d6', Bio = '#eb6834'), name = NULL) +
  labs(x = 'Method (ordered as in a)', y = 'Score (4 tissues)', tag = 'c',
       caption = 'Dashed: unintegrated data (median over tissues), blue batch, orange bio.') +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 8.5), legend.position = 'top',
        plot.tag = element_text(face = 'bold', size = 16), plot.caption = element_text(hjust = 0, colour = 'grey35'))

p10 = p10a / (p10b | p10c) + plot_layout(heights = c(1.25, 1)) +
  plot_annotation(caption = paste('Dissected datasets of four tissues. a: fill = rank among the 23 methods (darker = better);',
                                  'circle size = metric score; "–" = metric not defined for graph-based output;',
                                  'marker left of the name = input (blue PCA / own, orange UINMF, green scVI latent).'),
                  theme = theme(plot.caption = element_text(hjust = 0, size = 10, colour = 'grey30')))
save(p10, '10_summary_table.png', 17, 19)

# ---- 11. as 10, but every score is the gain over the raw unintegrated data of the same tissue ----
# dissected datasets only; per tissue, each metric (and Overall / Batch / Bio) minus that of the raw
# unintegrated data (PCA of shared genes), then averaged over the four tissues.
# Circle size = |gain|; red outline = below the unintegrated data. Fill = rank, as in 10.
raw_metric = expand_grid(tissue = names(tissue_labels)) |>
  pmap_dfr(function(tissue) list.files(score_dir(tissue, 'dissected'), pattern = '^05\\..*\\.csv$', full.names = TRUE) |>
             map_dfr(read_metrics) |> mutate(tissue = tissue)) |>
  filter(method == 'Unintegrated') |> select(tissue, metric, raw = value)
g_metric = diss_raw |> filter(metric %in% metric_cols$key) |>
  inner_join(raw_metric, by = c('tissue', 'metric')) |> mutate(value = value - raw) |>
  group_by(method, metric) |> summarise(value = mean(value, na.rm = TRUE), .groups = 'drop') |>
  mutate(value = if_else(is.nan(value), NA_real_, value))
raw_sum = scores |> filter(case == 'dissected', method == 'Unintegrated') |>
  select(tissue, overall0 = overall, batch0 = batch, bio0 = bio)
g_tissue = methods |> filter(case == 'dissected') |> left_join(raw_sum, by = 'tissue') |>
  mutate(d_overall = overall - overall0, d_batch = batch - batch0, d_bio = bio - bio0)
stopifnot(!anyNA(g_tissue$overall0))
g_sum = g_tissue |> group_by(method, input) |>
  summarise(Overall = mean(d_overall), Batch = mean(d_batch), Bio = mean(d_bio), .groups = 'drop') |>
  arrange(desc(Overall)) |> mutate(y = row_number())

# signed bars: one scale for the three bar columns, each with its zero placed to fit its range
g_bars = g_sum |> pivot_longer(c(Overall, Batch, Bio), names_to = 'key', values_to = 'value') |>
  left_join(layout_x, by = 'key') |> group_by(key) |> mutate(fill = rank_fill(value, group[1])) |>
  mutate(lo = min(0, value), hi = max(0, value)) |> ungroup()
k = 1.45 / max(g_bars |> group_by(key) |> summarise(r = max(hi) - min(lo)) |> pull(r))
g_bars = g_bars |> group_by(key) |> mutate(x0 = x + 0.05 - k * min(lo)) |> ungroup() |>
  mutate(xe = x0 + k * value)
g_zero = g_bars |> distinct(key, x0)
g_dots = g_metric |> rename(key = metric) |> inner_join(g_sum |> select(method, y), by = 'method') |>
  left_join(layout_x, by = 'key') |> group_by(key) |> mutate(fill = rank_fill(value, group[1])) |> ungroup() |>
  mutate(size = abs(value), loss = !is.na(value) & value < 0)
g_heads = layout_x |> mutate(xc = if_else(kind == 'bar', x + 0.75, x))
g_mark = g_sum |> mutate(col = input_cols[as.character(input)], shp = input_shapes[as.character(input)])
max_dot = max(abs(g_dots$value), na.rm = TRUE)

p11a = ggplot() +
  geom_segment(data = groups10, aes(x = x0, xend = x1, y = -0.55, yend = -0.55), linewidth = 0.5) +
  geom_text(data = groups10, aes(x = (x0 + x1) / 2, y = -1.35, label = label), fontface = 'bold', size = 4.2) +
  geom_text(data = g_heads, aes(x = xc, y = 0, label = label), size = 3.1, colour = '#52514e') +
  geom_segment(data = g_zero, aes(x = x0, xend = x0, y = 0.5, yend = nrow(g_sum) + 0.5), colour = 'grey55', linewidth = 0.4) +
  geom_rect(data = g_bars, aes(xmin = pmin(x0, xe), xmax = pmax(x0, xe), ymin = y - 0.33, ymax = y + 0.33, fill = fill),
            colour = 'grey35', linewidth = 0.25) +
  geom_text(data = g_bars, aes(x = pmax(x0, xe) + 0.06, y = y, label = sprintf('%+.2f', value)),
            hjust = 0, size = 2.5, colour = '#52514e') +
  geom_point(data = g_dots |> filter(!is.na(value)),
             aes(x = x, y = y, size = size, fill = fill, colour = loss, stroke = if_else(loss, 0.9, 0.3)), shape = 21) +
  geom_text(data = g_dots |> filter(is.na(value)), aes(x = x, y = y), label = '–', colour = 'grey60', size = 3) +
  geom_point(data = g_mark, aes(x = -0.45, y = y, shape = shp), colour = g_mark$col, size = 2.4) +
  geom_tile(data = leg, aes(gx, gy, fill = fill), width = 0.36, height = 0.75, colour = 'white') +
  annotate('text', x = 13.4 + 0.76, y = 0.9, label = 'Rank', fontface = 'bold', size = 3.3) +
  annotate('text', x = 14.85, y = 2.35, label = '1 (best)', hjust = 0, size = 2.8) +
  annotate('text', x = 14.85, y = 5.35, label = 'last', hjust = 0, size = 2.8) +
  scale_fill_identity() + scale_shape_identity() +
  scale_colour_manual(values = c(`FALSE` = 'grey30', `TRUE` = '#d03b3b'), guide = 'none') +
  scale_size_area(max_size = 7, limits = c(0, max_dot), name = '|Gain|') +
  scale_y_reverse(breaks = g_sum$y, labels = g_sum$method, expand = expansion(add = c(0.7, 0.3))) +
  scale_x_continuous(limits = c(-0.6, 15.6), expand = c(0, 0)) +
  guides(size = guide_legend(override.aes = list(fill = 'grey80'), title.position = 'top')) +
  coord_cartesian(clip = 'off') +
  labs(x = NULL, y = NULL, tag = 'a') +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = '#0b0b0b'),
        legend.position = c(0.965, 0.4), plot.tag = element_text(face = 'bold', size = 16),
        plot.margin = margin(30, 10, 5, 5))

# b: gain in bio vs gain in batch, mean +/- SE over the four tissues
g_sc = g_tissue |> mutate(bcol = case_when(str_starts(method, 'LIGER') ~ 'LIGER', TRUE ~ base)) |>
  group_by(method, input, bcol) |>
  summarise(mb = mean(d_batch), sb = sd(d_batch) / sqrt(n()), mo = mean(d_bio), so = sd(d_bio) / sqrt(n()), .groups = 'drop')
p11b = ggplot(g_sc, aes(mb, mo, colour = bcol, shape = input)) +
  geom_vline(xintercept = 0, colour = 'grey55') + geom_hline(yintercept = 0, colour = 'grey55') +
  geom_errorbar(aes(ymin = mo - so, ymax = mo + so), width = 0, alpha = 0.6) +
  geom_errorbar(aes(xmin = mb - sb, xmax = mb + sb), width = 0, alpha = 0.6, orientation = 'y') +
  geom_point(size = 3.2, stroke = 1.2) +
  scale_colour_manual(values = base_cols9, name = 'Integration method') +
  scale_shape_manual(values = input_shapes9, name = 'Input') +
  labs(x = 'Gain in batch correction (mean over 4 tissues)', y = 'Gain in bio conservation (mean over 4 tissues)',
       tag = 'b', caption = 'Bars: +/- 1 SE over tissues. Lines: no change from the unintegrated data.') +
  theme_classic(base_size = 11) +
  theme(legend.position = 'bottom', legend.box = 'vertical', plot.tag = element_text(face = 'bold', size = 16),
        plot.caption = element_text(hjust = 0, colour = 'grey35')) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(shape = 16)))

# c: gain in batch and bio per method over the four tissues
g_bx = g_tissue |> select(method, tissue, Batch = d_batch, Bio = d_bio) |>
  pivot_longer(c(Batch, Bio), names_to = 'component', values_to = 'gain') |>
  mutate(method = factor(method, levels = g_sum$method))
p11c = ggplot(g_bx, aes(method, gain, fill = component, colour = component)) +
  geom_hline(yintercept = 0, colour = 'grey55') +
  geom_boxplot(alpha = 0.2, width = 0.7, outlier.size = 0.8, linewidth = 0.4, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c(Batch = '#2a78d6', Bio = '#eb6834'), name = NULL) +
  scale_colour_manual(values = c(Batch = '#2a78d6', Bio = '#eb6834'), name = NULL) +
  labs(x = 'Method (ordered as in a)', y = 'Gain over unintegrated (4 tissues)', tag = 'c',
       caption = 'Line: no change from the unintegrated data.') +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 8.5), legend.position = 'top',
        plot.tag = element_text(face = 'bold', size = 16), plot.caption = element_text(hjust = 0, colour = 'grey35'))

# each panel in its own file (tags dropped)
no_tag = labs(tag = NULL)
save(p11a + no_tag + labs(caption = paste('Dissected datasets of four tissues; every score is minus that of the raw unintegrated',
                                          'data (PCA of shared genes) of the same tissue, then averaged over tissues.',
                                          '\nFill = rank among the 23 methods (darker = better); circle size = |gain|, red outline =',
                                          'below unintegrated; "–" = metric not defined for graph-based output; marker left of the name = input.')) +
       theme(plot.caption = element_text(hjust = 0, size = 10, colour = 'grey30')),
     '11a_summary_table_gain.png', 17, 11)
save(p11b + no_tag, '11b_gain_batch_vs_bio.png', 8.5, 8)
save(p11c + no_tag, '11c_gain_boxplots.png', 11, 6.5)
write_csv(g_tissue |> select(tissue, method, input, d_overall, d_batch, d_bio), file.path(out, 'gain_dissected.csv'))
print(g_sum |> mutate(across(c(Overall, Batch, Bio), \(x) round(x, 3))), n = Inf)
