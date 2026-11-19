# Appendix: bone marrow (dissected) UMAPs of five labellings of the same cells, 2 x 5 grid.
# Columns: unintegrated (joint PCA), scANVI, SAMap, BBKNN (native), Harmony on scVI-latent input; top row coloured by
# species, bottom row by cell type (EBCC labels, as in the benchmark). UMAPs from analysis/bm_umap_compute.py (identical settings,
# see figures/benchmark_summary/bm_umap/settings.txt). Same evaluation cells in every panel (those present in all five
# outputs), all drawn, in one random order shared by the panels. Species and cell-type colours as unintegrated_umap.png.
# No titles or footnotes in the image; caption in bm_umap_caption.txt.
# Usage: R_LIBS_SITE= R_LIBS_USER= <r45>/bin/Rscript --vanilla analysis/bm_umap_figure.R -> figures/benchmark_summary/bm_umap.{png,pdf}
setwd(Sys.getenv('PROJECT_ROOT', '.'))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(patchwork) })

# ---- config ----
umap_dir = 'figures/benchmark_summary/bm_umap'
labels_csv = 'figures/benchmark_summary/unintegrated_umap_data/bone_marrow.csv'   # cell, species, celltype (analysis/appendix_unintegrated_umap.R)
eval_file = 'bone_marrow/data/04_integrated_data/dissected_mouse/eval_cells.txt'
methods = c('Unintegrated' = 'Unintegrated', 'scANVI' = 'scANVI', 'SAMap' = 'SAMap', 'BBKNN' = 'BBKNN',
            'Harmony (scVI-latent)' = 'Harmony_scVI-latent')                         # column label = csv name
species_cols = c(human = '#4477AA', mouse = '#EE6677')                               # as analysis/visualisations.R
ct_palette = c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#E69F00', '#332288', '#8C564B',
               '#BBBBBB', '#44AA99', '#000000', '#9467BD', '#D62728', '#882255')     # as analysis/appendix_unintegrated_umap.R
zoom = c(0.01, 0.99)
point_size = 0.55
out = 'figures/benchmark_summary/bm_umap'
score_file = 'figures/benchmark_summary/benchmark_improvement_MAP_SAS.csv'   # scores shown under the column names
score_names = c('Unintegrated' = 'Unintegrated', 'scANVI' = 'scANVI', 'SAMap' = 'SAMap', 'BBKNN' = 'BBKNN',
                'Harmony (scVI-latent)' = 'Harmony scVI-latent')
graph_based = c('SAMap', 'BBKNN')

bm = read_csv(score_file, show_col_types = FALSE) |> filter(tissue == 'bone_marrow')
scores = bind_rows(bm |> filter(method %in% score_names) |> select(method, bio, batch, overall),
                   bm |> distinct(bio = ref_bio, batch = ref_batch, overall = ref_overall) |> mutate(method = 'Unintegrated'))
scores = scores[match(score_names, scores$method), ]
stopifnot(nrow(scores) == length(methods), !anyNA(scores$overall))
col_title = setNames(paste0(names(methods), ifelse(names(methods) %in% graph_based, ' \u2021', '')), names(methods))
col_sub = setNames(sprintf('bio %.3f · batch %.3f · overall %.3f', scores$bio, scores$batch, scores$overall), names(methods))
lab = read_csv(labels_csv, show_col_types = FALSE) |> select(cell, species, celltype) |>
  mutate(species = recode(species, hsapiens = 'human', mmusculus = 'mouse'))
eval_cells = readLines(eval_file)
# SAMap drops the species prefix; barcodes do not overlap between species, so the bare barcode identifies the cell
bare = sub('^(hsapiens|mmusculus)_', '', eval_cells)
stopifnot(!anyDuplicated(bare))
emb = lapply(methods, function(m) {
  e = read_csv(file.path(umap_dir, paste0(m, '.csv')), show_col_types = FALSE)
  i = match(e$cell, eval_cells); j = match(e$cell, bare)
  e$cell = eval_cells[ifelse(is.na(i), j, i)]
  e[!is.na(e$cell), ]
})
common = Reduce(intersect, lapply(emb, `[[`, 'cell'))
cat('evaluation cells:', length(eval_cells), '; in all five outputs:', length(common), '\n')
for (m in names(emb)) cat(sprintf('  %-22s %d evaluation cells\n', m, nrow(emb[[m]])))
set.seed(1)
order_cells = sample(common)                                                          # one draw order for every panel
types = sort(unique(lab$celltype[lab$cell %in% common]))
ct_cols = setNames(ct_palette[seq_along(types)], types)

d = bind_rows(lapply(names(emb), \(m) emb[[m]] |> mutate(method = m))) |>
  filter(cell %in% common) |> inner_join(lab, by = 'cell') |>
  mutate(method = factor(method, levels = names(methods)), ord = match(cell, order_cells)) |> arrange(ord)
stopifnot(all(table(d$method) == length(common)))

panel = function(m, col, values, name) {
  x = filter(d, method == m)
  xl = quantile(x$umap_1, zoom); yl = quantile(x$umap_2, zoom)
  xl = xl + c(-0.04, 0.04) * diff(xl); yl = yl + c(-0.04, 0.04) * diff(yl)
  ggplot(x, aes(umap_1, umap_2, colour = .data[[col]])) +
    ggrastr::rasterise(geom_point(size = point_size, stroke = 0), dpi = 300) +
    scale_colour_manual(values = values, name = name, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 3.5), nrow = if (col == 'species') 1 else 2)) +
    coord_equal(xlim = xl, ylim = yl) +
    labs(title = if (col == 'species') col_title[[m]] else NULL, subtitle = if (col == 'species') col_sub[[m]] else NULL) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = 'bold', size = 13, hjust = 0.5, margin = margin(b = 3)),
          plot.subtitle = element_text(size = 11, hjust = 0.5, colour = 'grey25', margin = margin(b = 6)),
          legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = 'bold'))
}
top = wrap_plots(lapply(names(methods), panel, col = 'species', values = species_cols, name = 'Species'), nrow = 1) +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
bottom = wrap_plots(lapply(names(methods), panel, col = 'celltype', values = ct_cols, name = 'Cell type'), nrow = 1) +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
p = wrap_elements(full = top) / wrap_elements(full = bottom)
ggsave(paste0(out, '.png'), p, width = 20, height = 10, dpi = 300, bg = 'white')
ggsave(paste0(out, '.pdf'), p, width = 20, height = 10, device = cairo_pdf, bg = 'white')
settings = readLines(file.path(umap_dir, 'settings.txt'))
writeLines(paste(sprintf(paste('Bone marrow (dissected mouse data), %s evaluation cells present in every output, all drawn in one random order.',
                               'Top: species; bottom: cell type (EBCC labels). Axes zoomed to the central 98%% of cells per panel.'),
                         format(length(common), big.mark = ',')),
                 paste('Under each name: the method\'s bio-conservation, batch-correction and overall scores in bone marrow (benchmark definition',
                       'with MAP and SAS; not improvements). \u2021 graph-based: no embedding, so bio = mean of NMI, ARI, cLISI and batch = mean of',
                       'iLISI, graph connectivity.'),
                 paste('UMAP settings, identical for all five panels and for the unintegrated UMAPs of the other appendix figure:', settings,
                       'A random initial layout is used throughout because spectral initialisation fails on some of these graphs.'), sep = '\n'), paste0(out, '_caption.txt'))
