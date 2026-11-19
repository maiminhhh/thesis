# Seurat label transfer with a same-species reference: how confident is it, and does its confidence mean anything?
# Run from the repository root. Reads the saved 02 annotation objects; nothing is refitted or rescored.
# Per tissue and species, all cells Seurat labels (celltype_seurat) with their prediction score (labelscore_seurat,
# the max prediction score from TransferData):
#   1. distribution of the score (bin width 0.05; a score of 1 in the last bin) and the share at exactly 1
#   2. cells split into fixed score bands; per band, the marker AUROC of Seurat's labels on those cells
#      (mean over the types Seurat calls in the band, marker score = mean expression of the type's markers)
# Usage:  Rscript analysis/seurat_confidence.R <tissue> <species>   -> figures/cross_tissue/seurat_confidence/<tissue>_<species>{,_hist}.csv
#         Rscript analysis/seurat_confidence.R                      figure -> figures/cross_tissue/seurat_confidence.{png,pdf,csv}

# helpers from analysis/vote_vs_ebcc.R (and through it analysis/visualisations.R)
ex = parse('analysis/vote_vs_ebcc.R')
stop_at = which(startsWith(sapply(ex, \(e) deparse(e)[1]), 'args = commandArgs'))
for (e in ex[seq_len(stop_at - 1)]) eval(e, globalenv())

out_dir = 'figures/cross_tissue/seurat_confidence'
bands = c(0, 0.5, 0.7, 0.9, 1 - 1e-9, Inf)
band_lab = c('< 0.5', '0.5–0.7', '0.7–0.9', '0.9–1', '= 1')

compute_one = function(tissue, species) {
  t0 = Sys.time()
  tick = function(msg) message(sprintf('[%s %s %5.1f min] %s', tissue, species, as.numeric(Sys.time() - t0, units = 'mins'), msg))
  a = tissues[[tissue]]$annotation[[species]]
  obj = readRDS(a$obj)
  assay = DefaultAssay(obj)
  rownames(obj) = obj[[assay]]@meta.data$gene_name
  X = LayerData(obj, assay = assay, layer = 'data')
  meta = obj@meta.data
  rm(obj); gc(); tick('object loaded')
  raw = as.character(meta[[annotators[['seurat']]]])
  lab = key(raw)
  s = meta$labelscore_seurat
  keep = which(!is.na(lab) & !is.na(s))
  markers = lapply(annotation_markers(a), intersect, rownames(X))
  markers = markers[lengths(markers) > 0]
  score = sapply(markers, \(g) Matrix::colMeans(X[g, keep, drop = FALSE]))   # kept cells x types
  colnames(score) = key(colnames(score))
  rm(X); gc()
  lab = lab[keep]; raw = raw[keep]; s = s[keep]
  band = cut(s, bands, labels = band_lab, right = FALSE)
  mean_auc_in = function(idx) {
    types = colnames(score)[vapply(colnames(score), \(k) { p = lab[idx] == k; any(p) && !all(p) }, logical(1))]
    if (!length(types)) return(c(NA_real_, 0))
    c(mean(vapply(types, \(k) auroc(score[idx, k], lab[idx] == k), numeric(1))), length(types))
  }
  res = purrr::map_dfr(band_lab, function(b) {
    idx = which(band == b)
    au = if (length(idx) >= 50) mean_auc_in(idx) else c(NA_real_, 0)   # too few cells for a stable AUROC
    tibble(band = b, n = length(idx), pct = length(idx) / length(s), auroc = au[1], n_types = au[2])
  }) |> mutate(tissue = tissue, species = species, n_cells = length(s), auroc_all = mean_auc_in(seq_along(s))[1],
               median_score = median(s), pct_one = mean(s >= 1 - 1e-9), .before = 1)
  # per cell type within each band (types with >= 20 of their cells in the band): the spread behind the mean
  band_types = purrr::map_dfr(band_lab, function(b) {
    idx = which(band == b)
    tt = intersect(colnames(score), names(which(table(lab[idx]) >= 20)))
    tt = tt[vapply(tt, \(k) !all(lab[idx] == k), logical(1))]
    tibble(band = b, type = tt, n_type = as.integer(table(lab[idx])[tt]),
           auroc = vapply(tt, \(k) auroc(score[idx, k], lab[idx] == k), numeric(1)))
  }) |> mutate(tissue = tissue, species = species, .before = 1)
  hist = tibble(bin_left = 0.05 * pmin(floor(s / 0.05 + 1e-9), 19)) |> dplyr::count(bin_left) |>
    mutate(tissue = tissue, species = species, pct = n / length(s), .before = 1)
  tick('done')
  # per cell type Seurat calls: score quantiles, share at exactly 1, and the type's marker AUROC
  # (its markers' score separating the cells Seurat gives the type from all other cells; NA without markers)
  types = tibble(type = raw, k = lab, s = s) |> group_by(type, k) |>
    summarise(n = n(), q10 = quantile(s, 0.1), q25 = quantile(s, 0.25), median = median(s),
              q75 = quantile(s, 0.75), q90 = quantile(s, 0.9), pct_one = mean(s >= 1 - 1e-9), .groups = 'drop') |>
    mutate(auroc = vapply(k, \(t) if (t %in% colnames(score)) auroc(score[, t], lab == t) else NA_real_, numeric(1)),
           tissue = tissue, species = species, pct_cells = n / length(s), .before = 1) |> select(-k)
  list(res = res, hist = hist, types = types, band_types = band_types)
}

args = commandArgs(TRUE)
if (length(args) == 2) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  r = compute_one(args[1], args[2])
  write_csv(r$res, file.path(out_dir, sprintf('%s_%s.csv', args[1], args[2])))
  write_csv(r$hist, file.path(out_dir, sprintf('%s_%s_hist.csv', args[1], args[2])))
  write_csv(r$types, file.path(out_dir, sprintf('%s_%s_types.csv', args[1], args[2])))
  write_csv(r$band_types, file.path(out_dir, sprintf('%s_%s_band_types.csv', args[1], args[2])))
  quit(save = 'no')
}

# ---- figure ----
cap_main = paste('Prediction score: Seurat TransferData\'s highest per-type score for the cell (0–1), from the same-species reference. All cells Seurat labels.',
                    'B: marker AUROC of Seurat\'s labels on the cells in each band (marker score = mean expression of the type\'s markers). Boxes: one AUROC per cell type with at least 20 cells',
                    'in the band (median, quartiles, whiskers to 1.5 IQR); diamond and line: mean over all the types Seurat calls in the band; dashed line: same mean on all cells.',
                 'Breast cancer mouse left blank: its reference has no Cancer Epithelial type, so its labels and scores are not comparable.', sep = '\n')
writeLines(cap_main, 'figures/cross_tissue/seurat_confidence_caption.txt')   # caption for the thesis text, not drawn
cap_type = paste('Per cell type Seurat assigns (same-species reference), all cells it labels. Open point: median prediction score; thick bar: 25–75%; thin bar: 10–90% of cells.',
                       'Right: cells assigned the type · % with a score of exactly 1 · marker AUROC of the type (its marker score separating the cells Seurat gives it from all others; – = no markers).',
                       'Types ordered by median score within each panel.',
                       'Breast cancer mouse left blank: its reference has no Cancer Epithelial type, so its labels and scores are not comparable.', sep = '\n')
writeLines(cap_type, 'figures/cross_tissue/seurat_confidence_by_type_caption.txt')   # caption for the thesis text, not drawn
# breast cancer mouse is left blank in both figures (panel kept, no data): reference lacks Cancer Epithelial
datasets = expand_grid(tissue = tissue_order, species = c('human', 'mouse')) |> filter(!(tissue == 'breast_cancer' & species == 'mouse'))
read_part = \(suffix) purrr::pmap_dfr(datasets, \(tissue, species)
  read_csv(file.path(out_dir, sprintf('%s_%s%s.csv', tissue, species, suffix)), show_col_types = FALSE))
tissue_name = \(t) str_to_sentence(str_replace(t, '_', ' '))
tissue_levels = tissue_name(tissue_order)
res = read_part('') |> mutate(tissue = factor(tissue_name(tissue), levels = tissue_levels), band = factor(band, levels = band_lab))
hist = read_part('_hist') |> mutate(tissue = factor(tissue_name(tissue), levels = tissue_levels))
write_csv(res, 'figures/cross_tissue/seurat_confidence.csv')
species_vals = c(human = species_cols[['hsapiens']], mouse = species_cols[['mmusculus']])

labs_one = distinct(res, tissue, species, pct_one, median_score, n_cells) |>
  mutate(label = sprintf('%s cells\nmedian %.2f\n%.0f%% at exactly 1', scales::comma(n_cells), median_score, 100 * pct_one))
pa = ggplot(hist, aes(bin_left + 0.025, pct, fill = species)) +
  geom_col(width = 0.05, alpha = 0.8, colour = 'white', linewidth = 0.2) +
  geom_text(data = labs_one, aes(x = 0.02, y = Inf, label = label), hjust = 0, vjust = 1.2, size = 2.8,
            colour = 'grey20', lineheight = 0.9, inherit.aes = FALSE) +
  facet_grid(species ~ tissue, labeller = labeller(species = str_to_sentence)) +
  scale_fill_manual(values = species_vals, guide = 'none') +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c('0', '0.25', '0.5', '0.75', '1')) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.45))) +
  labs(x = 'Seurat prediction score (bin width 0.05; the last bin includes 1)', y = 'Cells') +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = 'bold'),
        plot.title = element_text(face = 'bold', size = 11))
band_types = read_part('_band_types') |>
  mutate(tissue = factor(tissue_name(tissue), levels = tissue_levels), band = factor(band, levels = band_lab))
write_csv(band_types, 'figures/cross_tissue/seurat_confidence_band_types.csv')
dodge = position_dodge(width = 0.75)
pb = ggplot(band_types, aes(band, auroc, colour = species)) +
  geom_hline(data = distinct(res, tissue, species, auroc_all), aes(yintercept = auroc_all, colour = species, linetype = 'All cells (mean over types)'),
             linewidth = 0.5) +
  geom_boxplot(aes(group = interaction(band, species)), position = dodge, width = 0.65, outlier.size = 0.9,
               linewidth = 0.4, fill = NA) +
  geom_line(data = filter(res, !is.na(auroc)), aes(group = species, x = band, y = auroc), position = dodge, linewidth = 0.6, alpha = 0.7) +
  geom_point(data = filter(res, !is.na(auroc)), aes(shape = 'Mean over cell types'), position = dodge, size = 2.2, fill = 'white', stroke = 0.9) +
  facet_grid(~ tissue) +
  scale_colour_manual(values = species_vals, name = NULL, labels = str_to_sentence) +
  scale_shape_manual(values = 23, name = NULL) +
  scale_linetype_manual(values = 'dashed', name = NULL) +
  scale_y_continuous(breaks = seq(0.2, 1, 0.1)) +
  labs(x = 'Prediction score band', y = 'Marker AUROC of Seurat labels') +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 35, hjust = 1))
p = patchwork::wrap_plots(pa, pb, ncol = 1, heights = c(1.15, 1)) +
  patchwork::plot_annotation(
    tag_levels = 'A') &
  theme(legend.position = 'top')
ggsave('figures/cross_tissue/seurat_confidence.png', p, width = 13, height = 9.5, dpi = 300)
ggsave('figures/cross_tissue/seurat_confidence.pdf', p, width = 13, height = 9.5, device = cairo_pdf)
print(as.data.frame(res |> select(tissue, species, band, n, pct, auroc)), digits = 3)

# ---- figure 2: by cell type (the cell type Seurat assigns) ----
types = read_part('_types') |>
  mutate(tissue_lab = tissue_name(tissue), panel = sprintf('%s, %s', tissue_lab, species),
         panel = factor(panel, levels = as.vector(t(outer(tissue_levels, c('human', 'mouse'), \(a, b) sprintf('%s, %s', a, b))))),
         y = paste(type, panel, sep = '___'),
         label = sprintf('%s · %.0f%% at 1 · AUROC %s', scales::comma(n), 100 * pct_one, ifelse(is.na(auroc), '–', sprintf('%.2f', auroc))))
types = mutate(types, y = factor(y, levels = y[order(median, q25)]))   # by median score within each panel (free y scales)
write_csv(types, 'figures/cross_tissue/seurat_confidence_by_type.csv')
pt = ggplot(types, aes(y = y, colour = species)) +
  geom_linerange(aes(xmin = q10, xmax = q90), linewidth = 0.5, alpha = 0.6) +
  geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 2.2, alpha = 0.6) +
  geom_point(aes(x = median), size = 2, shape = 21, fill = 'white', stroke = 1) +
  geom_text(aes(x = 1.04, label = label), hjust = 0, size = 2.5, colour = 'grey25') +
  facet_wrap(~ panel, ncol = 2, scales = 'free_y', dir = 'h', drop = FALSE) +
  scale_y_discrete(labels = \(x) sub('___.*$', '', x)) +
  scale_colour_manual(values = species_vals, guide = 'none') +
  scale_x_continuous(limits = c(0, 1.72), breaks = seq(0, 1, 0.25), labels = c('0', '0.25', '0.5', '0.75', '1'),
                     expand = expansion(add = c(0.01, 0))) +
  coord_cartesian(clip = 'off') +
  labs(x = 'Seurat prediction score', y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_blank(),
        strip.text = element_text(face = 'bold', hjust = 0))
ggsave('figures/cross_tissue/seurat_confidence_by_type.png', pt, width = 13, height = 14, dpi = 300)
ggsave('figures/cross_tissue/seurat_confidence_by_type.pdf', pt, width = 13, height = 14, device = cairo_pdf)
cat('done figure\n')
