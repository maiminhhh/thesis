# Fetal liver, dissected dataset: are the leading methods for bio-conservation a group rather than a single winner?
# Leading methods = the union of the top 8 by bio-conservation under each labelling, plus any subsample winner (EBCC labels; best-labeller
# labels = human Seurat, mouse EBCC, the label-sensitivity rescoring), each scored under both labellings.
# Bio = mean of the six bio metrics (NMI, ARI, label ASW, isolated-label F1 and silhouette, cLISI), as the benchmark.
# Resampling: the same 50 independent 20,000-cell subsamples scored under each labelling (paired):
#   EBCC labels      bootstrap/fetal_liver/scores/<rep>/ (+ bootstrap_labels/ebcc/ for the two unintegrated inputs)
#   best labeller    bootstrap_labels/best_annotator/scores/<rep>/
# A: dot = score on all cells; interval = that score + (2.5% and 97.5% quantiles - median) of the method's
#    subsample scores under the same labelling. B: share of subsamples in which each method ranks first on bio
#    among all methods scored in that subsample, per labelling. Reads saved scores only.
# Usage: Rscript analysis/fetal_liver_leaders.R  -> figures/cross_tissue/fetal_liver_leaders.{png,pdf,csv}, _rank1.csv
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(readr); library(stringr); library(ggplot2) })
setwd(Sys.getenv('PROJECT_ROOT', '.'))

bio_metrics = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'isolated_label_F1', 'isolated_label_silhouette', 'cLISI')
n_top = 8
read_bio = function(files) {
  map_dfr(files, function(f) {
    df = read_csv(f, show_col_types = FALSE, name_repair = 'minimal')
    tibble(metric = df[[1]], value = suppressWarnings(as.numeric(df[[2]])), method = names(df)[2])
  }) |> filter(metric %in% bio_metrics, !is.na(value)) |>
    group_by(method) |> summarise(bio = mean(value), .groups = 'drop')   # mean of the available metrics
}
lab = c(ebcc = 'EBCC labels', alt = 'Best-labeller labels (human Seurat, mouse EBCC)')

# ---- full-data scores ----
dir = 'fetal_liver/data/05_evaluated_score/dissected_mouse'
full = bind_rows(
  mutate(read_bio(list.files(dir, '^05\\..*\\.csv$', full.names = TRUE)), labels = 'ebcc'),
  mutate(read_bio(list.files(file.path(dir, 'label_best_annotator'), '^05\\..*\\.csv$', full.names = TRUE)), labels = 'alt'))

# ---- subsample scores, per repeat and labelling ----
rep_files = \(dirs) map_dfr(dirs, \(d) tibble(f = list.files(d, '^05\\..*\\.csv$', full.names = TRUE, recursive = TRUE))) |>
  mutate(rep = as.integer(basename(dirname(f))))
boot = bind_rows(
  mutate(rep_files(c('bootstrap/fetal_liver/scores', 'bootstrap_labels/ebcc/scores')), labels = 'ebcc'),
  mutate(rep_files('bootstrap_labels/best_annotator/scores'), labels = 'alt')) |>
  group_by(labels, rep) |> group_modify(\(g, k) read_bio(g$f)) |> ungroup()

# leading methods: top 8 under either labelling (full data)
top = full |> group_by(labels) |> slice_max(bio, n = n_top) |> pull(method) |> unique()
# rank 1 is counted only on subsamples where every method is scored under both labellings (paired); the gaps are
# BBKNN-family scores lost to a transient scib graph-LISI read error, so the complete subsamples are a random subset
n_methods = n_distinct(boot$method)
ok_reps = boot |> count(rep) |> filter(n == 2 * n_methods) |> pull(rep)
# plus any method that ranks first in one of those subsamples (under either labelling)
winners = boot |> filter(rep %in% ok_reps) |> group_by(labels, rep) |> slice_max(bio, n = 1, with_ties = FALSE) |> pull(method) |> unique()
top = union(top, winners)
spread = boot |> filter(method %in% top) |> group_by(labels, method) |>
  summarise(n_rep = n(), d_lo = quantile(bio, 0.025) - median(bio), d_hi = quantile(bio, 0.975) - median(bio), .groups = 'drop')
d = full |> filter(method %in% top) |> inner_join(spread, by = c('labels', 'method')) |>
  mutate(lo = bio + d_lo, hi = bio + d_hi)
order = d |> filter(labels == 'ebcc') |> arrange(bio) |> pull(method)

# rank 1 in each complete subsample, among all methods
rank1 = boot |> filter(rep %in% ok_reps) |> group_by(labels, rep) |> slice_max(bio, n = 1, with_ties = FALSE) |> ungroup() |>
  count(labels, method, name = 'n_first') |>
  right_join(expand_grid(labels = names(lab), method = top), by = c('labels', 'method')) |>
  mutate(n_first = replace_na(n_first, 0L), share = n_first / length(ok_reps))

d = mutate(d, labels = factor(lab[labels], levels = lab), method = factor(method, levels = order))
rank1 = mutate(rank1, labels = factor(lab[labels], levels = lab), method = factor(method, levels = order))
write_csv(d, 'figures/cross_tissue/fetal_liver_leaders.csv')
write_csv(rank1, 'figures/cross_tissue/fetal_liver_leaders_rank1.csv')
print(as.data.frame(d |> select(method, labels, bio, lo, hi, n_rep) |> arrange(labels, desc(bio))), digits = 3)
print(as.data.frame(rank1 |> arrange(labels, desc(share))))
cat('repeats used for rank 1:', length(ok_reps), '\n')

cols = setNames(c('#2A7F9E', '#D2691E'), lab)
dodge = position_dodge(width = 0.55)
theme_fl = theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = 'bold', size = 11))
pa = ggplot(d, aes(bio, method, colour = labels, group = labels)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = 'y', linewidth = 0.8, position = dodge) +
  geom_point(size = 2.6, position = dodge) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = 'Bio-conservation score', y = NULL, title = 'Score on all cells, with 95% subsample interval') +
  theme_fl
pb = ggplot(rank1, aes(share, method, fill = labels)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = ifelse(n_first > 0, sprintf('%d/%d', n_first, length(ok_reps)), '')),
            position = position_dodge(width = 0.7), hjust = -0.15, size = 2.8, colour = 'grey25') +
  scale_fill_manual(values = cols, guide = 'none') +   # same colours as A's legend
  scale_x_continuous(labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = 'Share of subsamples ranked first', y = NULL, title = 'Ranked first on bio-conservation') +
  theme_fl + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p = patchwork::wrap_plots(pa, pb, widths = c(1.6, 1), guides = 'collect') +
  patchwork::plot_annotation(
    tag_levels = 'A',
    title = sprintf('Fetal liver, dissected mouse: leading methods by bio-conservation (top %d under either labelling, plus any subsample winner)', n_top),
    caption = paste('A: 95% intervals from 50 subsamples of 20,000 cells, the same cells scored under each labelling, centred on the full-data score.',
                    sprintf('B: share of the %d subsamples in which every method was scored under both labellings (BBKNN-family scores are missing in the rest) where the method ranks first among all methods.', length(ok_reps)),
                    'Each subsample holds 7% of the 282,334 cells, so the intervals overstate the uncertainty of the full-data scores. Methods ordered by the EBCC-label score.',
                    'Bio: mean of NMI, ARI, label ASW, isolated-label F1 and silhouette, cLISI (as the benchmark).', sep = '\n'),
    theme = theme(plot.title = element_text(face = 'bold', size = 12), plot.caption = element_text(hjust = 0, size = 8.5))) &
  theme(legend.position = 'top')
ggsave('figures/cross_tissue/fetal_liver_leaders.png', p, width = 12, height = 6, dpi = 300)
ggsave('figures/cross_tissue/fetal_liver_leaders.pdf', p, width = 12, height = 6, device = cairo_pdf)
