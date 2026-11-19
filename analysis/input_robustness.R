# Is the gain from the UINMF / scVI-latent input robust to which cells are scored?  Rscript analysis/input_robustness.R
# For every case (4 tissues x dissected + genes2/4/6/8) and every paired bootstrap subsample (20,000 cells):
# matched-algorithm differences, e.g. Harmony scVI-latent - Harmony, on the overall / bio / batch scores.
# Scores as analysis/benchmark_figures.R: bio = mean(NMI, ARI, ASW_label, cLISI [+ MAP]), batch = mean(PCR_batch,
# ASW_label/batch, iLISI, graph_conn [+ SAS]), overall = 0.6 bio + 0.4 batch; equal weights within each group.
# Uses the subsamples complete for every method of a case (works on a partial run).
# Bootstrap: dissected = bootstrap/summary/1_long.csv + bootstrap_extra/extra_long.csv; reduced genes =
# bootstrap_genes/genes_long.csv. Graph-only methods are not bootstrapped (their graph cannot be subsampled).
# Output: figures/input_robustness/*.{pdf,png,csv}; printout on stdout.
source('library.R')

# ---------------- config ----------------
out           = 'figures/input_robustness'
extra_bio     = c('MAP')       # as analysis/benchmark_figures.R; character(0) for the reference definition
extra_batch   = c('SAS')
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI')
algorithms    = c('Harmony', 'FastMNN', 'RPCA', 'CCA', 'Scanorama', 'LIGER')   # native method names
inputs        = c('scVI-latent', 'UINMF')
cases         = c(dissected = 0, genes2 = 20, genes4 = 40, genes6 = 60, genes8 = 80)   # % orthologues removed
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
level         = 0.95           # interval / decision level
min_rep       = 20             # fewer complete subsamples in a case: verdict 'pending'
top_k         = 5
# -----------------------------------------
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suffix = if (length(c(extra_bio, extra_batch))) paste0('_', paste(c(extra_bio, extra_batch), collapse = '_')) else ''

# ---- per-subsample metrics, all cases ----
dis = read_csv('bootstrap/summary/1_long.csv', show_col_types = FALSE) |>
  filter(`repeat` != 'point') |>
  transmute(tissue = dataset, case = 'dissected', rep = as.integer(`repeat`), method, metric, value = score,
            ref = score - delta_from_unintegrated)
dis_ref = dis |> group_by(tissue, case, rep, metric) |> summarise(value = mean(ref, na.rm = TRUE), .groups = 'drop') |>
  mutate(method = 'Unintegrated')
dis_extra = read_csv('bootstrap_extra/extra_long.csv', show_col_types = FALSE) |>
  pivot_longer(c(MAP, NC, SAS), names_to = 'metric', values_to = 'value') |>
  transmute(tissue, case = 'dissected', rep, method, metric, value)
gen = if (file.exists('bootstrap_genes/genes_long.csv'))
  read_csv('bootstrap_genes/genes_long.csv', show_col_types = FALSE) |>
    transmute(tissue, case = dataset, rep, method, metric, value = score) else NULL
met = bind_rows(select(dis, -ref), dis_ref, dis_extra, gen) |>
  filter(metric %in% c(batch_metrics, bio_metrics, extra_bio, extra_batch), !is.na(value),
         !str_detect(method, 'BBKNN|SAMap|scIB|^Unintegrated '))   # graph-only (not bootstrapped validly); unintegrated UINMF / scVI-latent inputs (no scib bootstrap)

sc = met |>
  mutate(group = ifelse(metric %in% c(batch_metrics, extra_batch), 'batch', 'bio')) |>
  group_by(tissue, case, rep, method, group) |> summarise(v = mean(value), n = n(), .groups = 'drop') |>
  pivot_wider(names_from = group, values_from = c(v, n)) |>
  rename(bio = v_bio, batch = v_batch) |>
  mutate(overall = 0.6 * bio + 0.4 * batch)
# a subsample counts only when every method of the case has all its metrics
need = c(length(bio_metrics) + length(extra_bio), length(batch_metrics) + length(extra_batch))
n_methods = sc |> group_by(tissue, case) |> summarise(m = n_distinct(method), .groups = 'drop')
complete = sc |> filter(n_bio == need[1], n_batch == need[2]) |> dplyr::count(tissue, case, rep, name = 'k') |>
  inner_join(n_methods, by = c('tissue', 'case')) |>
  group_by(tissue, case) |> filter(k == max(m)) |> ungroup() |> select(tissue, case, rep)
sc = semi_join(sc, complete, by = c('tissue', 'case', 'rep'))
reps = dplyr::count(complete, tissue, case, name = 'n_rep')
cat('Complete subsamples per case:\n'); print(as.data.frame(pivot_wider(reps, names_from = case, values_from = n_rep)), row.names = FALSE)

# ---- matched-algorithm differences ----
pairs = expand_grid(algorithm = algorithms, input = inputs) |> mutate(method = paste(algorithm, input)) |>
  filter(method %in% unique(sc$method))
long = sc |> pivot_longer(c(overall, bio, batch), names_to = 'score', values_to = 'value') |>
  select(tissue, case, rep, method, score, value)
d = pairs |> inner_join(long, by = 'method') |>
  inner_join(long |> rename(algorithm = method, native = value), by = c('algorithm', 'tissue', 'case', 'rep', 'score')) |>
  mutate(diff = value - native)
a = (1 - level) / 2
summ = d |> group_by(algorithm, input, tissue, case, score) |>
  summarise(n_rep = n(), median = median(diff), lo = quantile(diff, a), hi = quantile(diff, 1 - a),
            p_better = mean(diff > 0), .groups = 'drop') |>
  mutate(verdict = case_when(n_rep < min_rep ~ 'pending', lo > 0 ~ 'better', hi < 0 ~ 'worse', TRUE ~ 'unresolved'),
         removed = cases[case])
write_csv(summ, file.path(out, paste0('matched_differences', suffix, '.csv')))

cat(sprintf('\n== Matched algorithm: input vs native, overall score; cases where the %d%% interval is above / below 0 ==\n',
            round(100 * level)))
summ |> filter(score == 'overall') |> group_by(input, algorithm) |>
  summarise(cases = n(), better = sum(verdict == 'better'), worse = sum(verdict == 'worse'),
            unresolved = sum(verdict == 'unresolved'), pending = sum(verdict == 'pending'),
            median_diff = round(median(median), 3), .groups = 'drop') |>
  arrange(input, -better) |> as.data.frame() |> print(row.names = FALSE)
cat('\nSame, by score and input (all algorithms):\n')
summ |> group_by(score, input) |>
  summarise(comparisons = n(), better = sum(verdict == 'better'), worse = sum(verdict == 'worse'),
            unresolved = sum(verdict == 'unresolved'), pending = sum(verdict == 'pending'), .groups = 'drop') |>
  as.data.frame() |> print(row.names = FALSE)

# ---- group advantage: mean over algorithms of (input - native), per subsample ----
grp = d |> group_by(input, tissue, case, rep, score) |> summarise(adv = mean(diff), .groups = 'drop') |>
  group_by(input, tissue, case, score) |>
  summarise(n_rep = n(), median = median(adv), lo = quantile(adv, a), hi = quantile(adv, 1 - a),
            p_better = mean(adv > 0), .groups = 'drop') |>
  mutate(removed = cases[case])
write_csv(grp, file.path(out, paste0('group_advantage', suffix, '.csv')))
cat('\n== Mean advantage over the same algorithm on native input (overall), per case: median [interval] ==\n')
grp |> filter(score == 'overall', n_rep >= min_rep) |>
  mutate(txt = sprintf('%+.3f [%+.3f, %+.3f]', median, lo, hi)) |>
  select(input, tissue, case, txt) |> pivot_wider(names_from = case, values_from = txt) |>
  as.data.frame() |> print(row.names = FALSE)

# ---- top-k membership per case ----
topk = sc |> filter(method != 'Unintegrated') |> group_by(tissue, case, rep) |>
  mutate(r = rank(-overall)) |> group_by(tissue, case, method) |>
  summarise(p_top = mean(r <= top_k), .groups = 'drop')
write_csv(topk, file.path(out, paste0('top', top_k, '_frequency', suffix, '.csv')))
cat(sprintf('\n== Share of subsamples in the top %d (overall), mean over cases, and cases with share >= 0.5 ==\n', top_k))
topk |> group_by(method) |> summarise(mean_share = round(mean(p_top), 2), cases_majority = sum(p_top >= 0.5), cases = n()) |>
  arrange(-mean_share) |> head(10) |> as.data.frame() |> print(row.names = FALSE)

# ---- figure 1: matched differences, every case ----
tl = function(x) factor(tissue_labels[x], levels = tissue_labels)
f1 = summ |> mutate(row = factor(paste0(algorithm, ' (', input, ')'),
                                 levels = rev(c(paste0(algorithms, ' (scVI-latent)'), paste0(algorithms, ' (UINMF)')))),
                    col = factor(paste0(removed, '%'), levels = paste0(cases, '%')), tissue = tl(tissue),
                    score = factor(score, levels = c('overall', 'bio', 'batch'),
                                   labels = c('Overall (0.6 bio + 0.4 batch)', 'Bio-conservation', 'Batch correction')),
                    mark = c(better = '+', worse = '−', unresolved = '', pending = '·')[verdict]) |>
  filter(!is.na(row))
lim = max(abs(f1$median))
p1 = ggplot(f1, aes(col, row)) +
  geom_tile(aes(fill = median), colour = 'white', linewidth = 0.5) +
  geom_text(aes(label = mark), size = 3.6, fontface = 'bold') +
  facet_grid(score ~ tissue) +
  scale_fill_gradient2(low = '#b2182b', mid = '#f7f7f7', high = '#2166ac', midpoint = 0, limits = c(-lim, lim),
                       name = 'Median difference\n(input - native)') +
  labs(x = 'Shared orthologues removed', y = NULL) +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = 'bold'), legend.position = 'right')
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('matched_differences', suffix, '.', ext)), p1, width = 12, height = 10, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')

# ---- figure 2: group advantage vs orthologue removal, with bootstrap interval ----
input_cols = c(`scVI-latent` = '#1baf7a', UINMF = '#eb6834')
f2 = grp |> mutate(tissue = tl(tissue), score = factor(score, levels = c('overall', 'bio', 'batch'),
                   labels = c('Overall', 'Bio-conservation', 'Batch correction')))
p2 = ggplot(f2, aes(removed, median, colour = input, fill = input)) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40', linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) + geom_point(size = 1.8) +
  facet_grid(score ~ tissue, scales = 'free_y') +
  scale_x_continuous(breaks = cases, labels = \(x) paste0(x, '%')) +
  scale_colour_manual(values = input_cols, name = 'Input') + scale_fill_manual(values = input_cols, name = 'Input') +
  labs(x = 'Shared orthologues removed (0% = dissected)', y = 'Advantage over the same algorithm on native input') +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = 'grey95', colour = 'grey70'),
        strip.text = element_text(face = 'bold'), legend.position = 'bottom')
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('group_advantage', suffix, '.', ext)), p2, width = 13, height = 8, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')
