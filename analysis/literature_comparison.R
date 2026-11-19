# Where do published scRNA-seq cross-species benchmarks place the basic methods, relative to the spread of our tissues?
#   Rscript analysis/literature_comparison.R
# Ours: overall improvement over unintegrated, scib metrics only (no MAP / SAS), dissected; basic methods = every
# algorithm on its native input + LIGER UINMF; ranked within each tissue.
# Published (human-mouse tasks only, one-to-one orthologue input):
#   Song et al. 2023 (BENGAL), Supplementary Data 1: Pancreas_hs_mm, Heart_hs_mm; O2O homology mapping; integrated_score.
#   Zhong et al. 2025, Supplementary Table S9: task4, task21 (human-mouse, cross order); Overall.Score.
#   Wang et al. 2026, Supplementary Table S6: task1 (human-mouse, O2O genes); Predicted.Overall.Score (their model).
# Rank percentile = (rank - 1) / (n - 1), 0 = best: ours among the basic methods in each tissue; a published task's
# among the methods it shares with our basic set. One point per tissue / per published task.
# Rows: ordered by the mean published rank percentile.
# Output: figures/literature/literature_comparison.{pdf,png,csv}; printout.
source('library.R')

# ---------------- config ----------------
ours_file  = 'figures/benchmark_summary/benchmark_improvement.csv'   # reference definition (scib only)
out        = 'figures/literature'
song_file  = 'scrna_benchmark/41467_2023_41855_MOESM4_ESM.csv'
zhong_file = 'scrna_benchmark/Supplementary Tables.xlsx'
wang_file  = 'scrna_benchmark/supplementary_tables_bbag490 (1).xlsx'
song_tasks = c('Pancreas_hs_mm', 'Heart_hs_mm');  song_homology = 'O2O'
zhong_tasks = c('task4', 'task21')
wang_tasks  = c('task1')
name_map = c(fastMNN = 'FastMNN', `SeuratV4 CCA` = 'CCA', `SeuratV4 RPCA` = 'RPCA', `Seurat v4 CCA` = 'CCA',
             `Seurat v5 CCA` = 'CCA', `Seurat v5 RPCA` = 'RPCA')
paper_labels = c(song = 'Song et al. 2023 (BENGAL)', zhong = 'Zhong et al. 2025', wang = 'Wang et al. 2026')
paper_shapes = c(song = 24, zhong = 22, wang = 23)
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
tissue_cols   = c(bone_marrow = '#E69F00', colon = '#009E73', breast_cancer = '#CC79A7', fetal_liver = '#56B4E9')   # as analysis/benchmark_figures.R
# -----------------------------------------
dir.create(out, showWarnings = FALSE, recursive = TRUE)
pct = function(score) { r = rank(-score); (r - 1) / (length(r) - 1) }
std = function(m) ifelse(m %in% names(name_map), name_map[m], m)

ours = read_csv(ours_file, show_col_types = FALSE) |> filter(input == 'Native' | method == 'LIGER UINMF') |>
  select(method, tissue, score = imp_overall)
basic = sort(unique(ours$method))
ours = ours |> group_by(tissue) |> mutate(rank = rank(-score), pct = pct(score)) |> ungroup()

song = read_csv(song_file, skip = 2, show_col_types = FALSE) |>
  filter(task_name %in% song_tasks, homology_method_batch == song_homology) |>
  transmute(paper = 'song', task = task_name, method = std(integration_method_batch), score = integrated_score)
zhong = readxl::read_excel(zhong_file, sheet = 'Supplementary Table S9', skip = 1) |>
  filter(Task %in% zhong_tasks) |> transmute(paper = 'zhong', task = Task, method = std(Method), score = as.numeric(Overall.Score))
wang = readxl::read_excel(wang_file, sheet = 'Supplementary Table S6', skip = 1) |>
  filter(Task %in% wang_tasks) |> transmute(paper = 'wang', task = Task, method = std(Method), score = as.numeric(Predicted.Overall.Score))
pub_all = bind_rows(song, zhong, wang)
cat('Published rows per task (all methods) and methods shared with our basic set:\n')
pub_all |> group_by(paper, task) |> summarise(n_all = n(), n_na = sum(is.na(score)),
                                              shared = paste(sort(intersect(method, basic)), collapse = ', '), .groups = 'drop') |>
  as.data.frame() |> print(row.names = FALSE)
pub = pub_all |> filter(method %in% basic, !is.na(score)) |> group_by(paper, task) |>
  mutate(n = n(), rank = rank(-score), pct = pct(score)) |> ungroup()

# ---- table ----
tab = ours |> group_by(method) |> summarise(mean_rank = mean(rank), tissue_min = min(pct), tissue_max = max(pct)) |>
  left_join(pub |> mutate(col = paste(paper, task, sep = ':')) |> select(method, col, pct) |>
              pivot_wider(names_from = col, values_from = pct), by = 'method') |> arrange(mean_rank)
write_csv(tab, file.path(out, 'literature_comparison.csv'))
pub_cols = setdiff(names(tab), c('method', 'mean_rank', 'tissue_min', 'tissue_max'))
tab = tab |> mutate(n_pub = rowSums(!is.na(across(all_of(pub_cols)))),
                    n_inside = rowSums(across(all_of(pub_cols), \(p) !is.na(p) & p >= tissue_min - 1e-9 & p <= tissue_max + 1e-9)))
cat('\nRank percentile (0 = best). Ours: range over tissues; published: per task among shared methods.\n')
tab |> mutate(across(where(is.double), \(x) round(x, 2))) |> as.data.frame() |> print(row.names = FALSE)
cat(sprintf('\nPublished positions inside the range of our tissues: %d of %d\n', sum(tab$n_inside), sum(tab$n_pub)))

# ---- figure ----
# rows ordered by the literature: mean published rank percentile over all human-mouse tasks (each task counts once)
lit = pub |> group_by(method) |> summarise(lit_mean = mean(pct), n_tasks = n()) |> arrange(lit_mean)
cat('\nRow order (mean published rank percentile):\n'); print(as.data.frame(lit |> mutate(lit_mean = round(lit_mean, 2))), row.names = FALSE)
ord = c(lit$method, setdiff(tab$method, lit$method))
lab = function(m) ifelse(m %in% c('CCA', 'RPCA'), paste0(m, ' (Seurat)'), m)
yl = setNames(rev(seq_along(ord)), ord)
o = ours |> mutate(y = yl[method] + 0.16, tissue = factor(tissue, levels = names(tissue_labels)))
rng = o |> group_by(method, y) |> summarise(lo = min(pct), hi = max(pct), .groups = 'drop')
p = pub |> mutate(paper = factor(paper, levels = names(paper_labels))) |>
  mutate(y = yl[method] - 0.2 + (as.integer(paper) - 2) * 0.09)
fig = ggplot() +
  geom_hline(yintercept = seq_along(ord), colour = 'grey93', linewidth = 4) +
  geom_segment(data = rng, aes(x = lo, xend = hi, y = y, yend = y), colour = 'grey55', linewidth = 0.6) +
  geom_point(data = o, aes(pct, y, colour = tissue), size = 2.4) +
  geom_point(data = p, aes(pct, y, shape = paper), size = 2.6, fill = 'grey25', colour = 'black', stroke = 0.4) +
  scale_y_continuous(breaks = yl, labels = lab(names(yl)), expand = expansion(add = 0.6)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c('0\nbest', '0.25', '0.5', '0.75', '1\nworst')) +
  scale_colour_manual(values = tissue_cols, labels = tissue_labels, name = 'This study (tissue)') +
  scale_shape_manual(values = paper_shapes, labels = paper_labels, name = 'scRNA-seq benchmark\n(human-mouse task)') +
  labs(x = 'Rank percentile', y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.y = element_text(colour = 'grey10'),
        legend.position = 'right')
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('literature_comparison.', ext)), fig, width = 9.5, height = 6.5, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')

# ---- heatmap: our tissues, our mean, one column per published human-mouse task ----
# Fill = rank percentile; number = rank (ours among the 11 basic methods; a task's among the methods its paper shares
# with us, printed as rank/n; .5 = tie). Rows: our mean rank over tissues.
task_cols = c(`song:Pancreas_hs_mm` = 'BENGAL\npancreas', `song:Heart_hs_mm` = 'BENGAL\nheart',
              `zhong:task4` = 'Zhong 4\nmulti-tissue', `zhong:task21` = 'Zhong 21\neye', `wang:task1` = 'Wang 1\nbrain')
# Zhong (Tables S2-S3): task 4 = adult multi-tissue atlases (human 5 tissues, GSE134355; mouse 6 tissues, GSE153562);
# task 21 = eye, aqueous humour outflow (SCP780), 11 shared cell types, 1,323 mouse cells.
# Wang (Tables S1, S4): task 1 = brain: human (GSE186538) + mouse (GSE158450, isoform atlas of the postnatal mouse
# brain), 10x, 39,630 cells.
n_basic = length(basic)
hm_ours = bind_rows(ours |> transmute(method, col = tissue_labels[tissue], pct, rank = as.character(rank)),
                    ours |> group_by(method) |> summarise(pct = mean(pct), rank = sprintf('%.1f', mean(rank))) |> mutate(col = 'Mean'))
ties = pub |> group_by(paper, task) |> summarise(tied = sum(rank %% 1 != 0), .groups = 'drop')
cat('\nTied ranks within each published task:\n'); print(as.data.frame(ties), row.names = FALSE)
hdr = unname(task_cols)
hm_pub = pub |> transmute(method, col = task_cols[paste(paper, task, sep = ':')], pct,
                          rank = paste0(ifelse(rank %% 1 == 0, sprintf('%d', as.integer(rank)), sprintf('%.1f', rank)), '/', n))
stopifnot(!anyNA(hm_pub$col))
for (pp in c('song', 'zhong')) {   # agreement between a paper's two human-mouse tasks
  w = pub |> filter(paper == pp) |> select(task, method, rank) |> pivot_wider(names_from = task, values_from = rank)
  cat(sprintf('Kendall tau between the two %s tasks: %.2f\n', pp, cor(w[[2]], w[[3]], method = 'kendall')))
}
cols = c(tissue_labels, Mean = 'Mean', hdr)
gap = c(rep(0, 4), 0.25, 0.5, 0.5, 0.75, 0.75, 1)   # gaps: before Mean, before the benchmarks, between papers
x_of = setNames(seq_along(cols) + gap, cols)
ord2 = tab$method   # heatmap rows: our mean rank over tissues, best at the top
yl2 = setNames(rev(seq_along(ord2)), ord2)
hm = bind_rows(hm_ours, hm_pub) |> complete(method = ord2, col = cols) |> mutate(x = x_of[col], y = yl2[method])
fig2 = ggplot(hm, aes(x, y)) +
  geom_tile(aes(fill = pct), colour = 'white', linewidth = 0.8, width = 0.96, height = 0.94) +
  annotate('rect', xmin = x_of[['Mean']] - 0.5, xmax = x_of[['Mean']] + 0.5, ymin = 0.5, ymax = max(yl2) + 0.5,
           fill = NA, colour = 'grey20', linewidth = 0.7) +
  geom_text(data = filter(hm, is.na(pct)), label = '–', colour = 'grey55', size = 3.6) +
  geom_text(data = filter(hm, !is.na(pct)), aes(label = rank, colour = pct < 0.45), size = 3.4,
            fontface = ifelse(filter(hm, !is.na(pct))$col == 'Mean', 'bold', 'plain')) +
  scale_fill_gradient(low = '#3f007d', high = '#f2f0f7', limits = c(0, 1), na.value = 'grey92', breaks = c(0, 0.5, 1),
                      labels = c('0 (best)', '0.5', '1 (worst)'), name = 'Rank percentile') +
  scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = 'grey10'), guide = 'none') +
  scale_x_continuous(breaks = unname(x_of), labels = names(x_of), position = 'top', expand = expansion(add = 0.1)) +
  scale_y_continuous(breaks = yl2, labels = lab(names(yl2)), expand = expansion(add = 0.1)) +
  annotate('text', x = mean(x_of[c(tissue_labels, 'Mean')]), y = max(yl2) + 2.15, label = sprintf('This study (rank of %d)', n_basic),
           fontface = 'bold', size = 3.6) +
  annotate('text', x = mean(x_of[hdr]), y = max(yl2) + 2.15, label = 'Human-mouse scRNA-seq benchmarks\n(rank among shared methods)',
           fontface = 'bold', size = 3.6, lineheight = 0.9) +
  coord_cartesian(clip = 'off', ylim = c(0.5, max(yl2) + 0.5)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = 'grey10'), axis.text.x.top = element_text(size = 9.5),
        legend.position = 'bottom', plot.margin = margin(60, 10, 5, 5))
for (ext in c('pdf', 'png'))
  ggsave(file.path(out, paste0('literature_heatmap.', ext)), fig2, width = 11.5, height = 6.5, dpi = 300, bg = 'white',
         device = if (ext == 'pdf') cairo_pdf else 'png')
