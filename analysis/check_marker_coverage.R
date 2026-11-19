# How many of the specified marker genes are actually measured, and expressed, in each tissue?
# Per tissue, species and cell type: markers specified in the annotation script (annotation_markers), those in the
# dataset's gene panel (post-QC object), and of those the ones detected (count > 0) in >= 10% of the cells Seurat
# (same-species reference) assigns the type. Run from the repository root; reads saved objects only.
# Usage: Rscript analysis/check_marker_coverage.R <tissue> <species> -> figures/cross_tissue/marker_coverage/<tissue>_<species>.csv
#        Rscript analysis/check_marker_coverage.R                    -> figures/cross_tissue/marker_coverage.csv, printed summary
ex = parse('analysis/vote_vs_ebcc.R')
stop_at = which(startsWith(sapply(ex, \(e) deparse(e)[1]), 'args = commandArgs'))
for (e in ex[seq_len(stop_at - 1)]) eval(e, globalenv())

out_dir = 'figures/cross_tissue/marker_coverage'
min_detect = 0.10
args = commandArgs(TRUE)
if (length(args) == 2) {
  dir.create(out_dir, showWarnings = FALSE)
  a = tissues[[args[1]]]$annotation[[args[2]]]
  obj = readRDS(a$obj)
  assay = DefaultAssay(obj)
  rownames(obj) = obj[[assay]]@meta.data$gene_name
  X = LayerData(obj, assay = assay, layer = 'counts')
  lab = key(as.character(obj@meta.data[[annotators[['seurat']]]]))
  rm(obj); gc()
  markers = annotation_markers(a)
  res = purrr::imap_dfr(markers, function(g, type) {
    g = unique(g); inp = intersect(g, rownames(X))
    cells = which(lab == key(type))
    det = if (length(cells) && length(inp)) Matrix::rowMeans(X[inp, cells, drop = FALSE] > 0) else setNames(numeric(0), character(0))
    tibble(type = type, n_specified = length(g), n_in_panel = length(inp),
           n_detected = sum(det >= min_detect), n_cells_seurat = length(cells),
           missing = paste(setdiff(g, rownames(X)), collapse = ';'),
           in_panel_not_detected = paste(names(det)[det < min_detect], collapse = ';'))
  }) |> mutate(tissue = args[1], species = args[2], n_panel_genes = nrow(X), .before = 1)
  write_csv(res, file.path(out_dir, sprintf('%s_%s.csv', args[1], args[2])))
  quit(save = 'no')
}
d = purrr::map_dfr(list.files(out_dir, '\\.csv$', full.names = TRUE), read_csv, show_col_types = FALSE) |>
  mutate(tissue = factor(tissue, levels = tissue_order)) |> arrange(tissue, species, type)
write_csv(d, 'figures/cross_tissue/marker_coverage.csv')
summ = d |> group_by(tissue, species) |>
  summarise(types = n(), specified = sum(n_specified), in_panel = sum(n_in_panel),
            pct_in_panel = 100 * in_panel / specified, detected = sum(n_detected),
            types_lt3_in_panel = sum(n_in_panel < 3), types_none_detected = sum(n_detected == 0 & n_cells_seurat > 0),
            .groups = 'drop')
print(as.data.frame(summ), digits = 3)
print(as.data.frame(d |> select(tissue, species, type, n_specified, n_in_panel, n_detected, n_cells_seurat)), width = 200)

# ---- figure: how many of each panel's genes are markers, per cell type ----
# panel genes = the raw panel's gene targets (figures/cross_tissue/raw_gene_detection.csv, analysis/check_zero_genes.py; same
# genes as the post-QC objects); markers = annotation_markers. Breast cancer mouse left blank (panel kept).
pretty_type = \(x) x |> str_replace('^Epithelial_', 'Epithelial: ') |> str_replace('Stem_Progenitor', 'Stem/Progenitor') |>
  str_replace('^B_Plasma$', 'B/Plasma') |> str_replace('^Neural_Glial$', 'Neural/Glial') |>
  str_replace('^Stromal_Mesenchymal$', 'Stromal/Mesenchymal') |> str_replace('^Smooth_Muscle$', 'Smooth muscle')
panel_genes = read_csv('figures/cross_tissue/raw_gene_detection.csv', show_col_types = FALSE) |> filter(feature_type == 'Gene Expression')
tissue_name = \(t) str_to_sentence(str_replace(t, '_', ' '))
counts = purrr::pmap_dfr(expand_grid(tissue = tissue_order, species = c('human', 'mouse')), function(tissue, species) {
  panel = panel_genes$feature[panel_genes$tissue == tissue & panel_genes$species == species]
  m = lapply(annotation_markers(tissues[[tissue]]$annotation[[species]]), \(g) intersect(unique(g), panel))
  tibble(tissue = tissue, species = species, type = pretty_type(names(m)), n_markers = lengths(m),
         n_panel = length(panel), n_any = length(unique(unlist(m))))
})
stopifnot(all(counts$n_markers == (d |> arrange(match(tissue, tissue_order)) |> pull(n_in_panel))[match(
  paste(counts$tissue, counts$species, counts$type), paste(d$tissue, d$species, pretty_type(d$type)))]))
write_csv(counts, 'figures/cross_tissue/panel_markers_per_type.csv')
heads = counts |> distinct(tissue, species, n_panel, n_any) |>
  mutate(panel = sprintf('%s, %s: %d of %d panel genes are markers', tissue_name(tissue), species, n_any, n_panel),
         panel = if_else(tissue == 'breast_cancer' & species == 'mouse', sprintf('%s, %s', tissue_name(tissue), species), panel))
fig = counts |> inner_join(select(heads, tissue, species, panel), by = c('tissue', 'species')) |>
  mutate(panel = factor(panel, levels = heads$panel), y = paste(type, panel, sep = '___'))
fig = mutate(fig, y = factor(y, levels = y[order(panel, n_markers, decreasing = c(TRUE, FALSE), method = 'radix')]))
shown = filter(fig, !(tissue == 'breast_cancer' & species == 'mouse'))
p = ggplot(shown, aes(n_markers, y)) +
  geom_col(width = 0.7, fill = '#2A7F9E') +
  geom_text(aes(label = n_markers), hjust = -0.3, size = 2.7, colour = 'grey25') +
  facet_wrap(~ panel, ncol = 2, scales = 'free', drop = FALSE) +
  scale_y_discrete(labels = \(x) sub('___.*$', '', x)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12)), breaks = scales::breaks_pretty(4)) +
  labs(x = 'Panel genes that are markers of the cell type', y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = 'bold', hjust = 0))
ggsave('figures/cross_tissue/marker_coverage.png', p, width = 12, height = 13, dpi = 300)
ggsave('figures/cross_tissue/marker_coverage.pdf', p, width = 12, height = 13, device = cairo_pdf)
writeLines(paste('Number of genes in each dataset\'s gene panel that are markers of each cell type in the annotation (used for marker AUROC).',
                 'Panel headings: panel genes that are a marker of at least one type (a gene can mark several types), out of all panel genes.',
                 'Breast cancer mouse left blank: its reference has no Cancer Epithelial type.', sep = '\n'),
           'figures/cross_tissue/marker_coverage_caption.txt')
print(as.data.frame(distinct(counts, tissue, species, n_panel, n_any)))
cat('done figure\n')
