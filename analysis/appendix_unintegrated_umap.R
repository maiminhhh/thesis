# Appendix: UMAP of the unintegrated data (joint log-normalised PCA, the reference of the benchmark), dissected mouse case,
# one row per tissue: species | cell type, labelled with the unintegrated bio / batch / overall scores
# (figures/benchmark_summary/benchmark_improvement<score_suffix>.csv, ref_* columns).
# UMAP recomputed from the saved joint PCA with the same settings as the bone-marrow method figure (bm_umap.png), by
# analysis/appendix_unintegrated_umap_compute.py -> unintegrated_umap_data/<tissue>.csv (run that first); evaluation cells only,
# species = orig.ident, cell type = celltype. (Before 2026-11-05 this used the Seurat RunUMAP X_umap saved in the h5ad.)
# Usage: R_LIBS_SITE= R_LIBS_USER= <r45>/bin/Rscript --vanilla analysis/appendix_unintegrated_umap.R
#   -> figures/benchmark_summary/unintegrated_umap.{png,pdf}, _scores.csv, _caption.txt (PDF: point layers rasterised with ggrastr)
setwd(Sys.getenv('PROJECT_ROOT', '.'))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(patchwork) })

# ---- config ----
cases = c(colon         = 'colon/data/03_concatenated_data/03_concatenated_data_subset|colon/data/04_integrated_data/dissected_mouse/',
          bone_marrow   = 'bone_marrow/data/03_concatenated_data/03_concatenated_data_subset|bone_marrow/data/04_integrated_data/dissected_mouse/',
          fetal_liver   = 'fetal_liver/data/03_concatenated_data/03_concatenated_data_subset|fetal_liver/data/04_integrated_data/dissected_mouse/',
          breast_cancer = 'breast_cancer/data/dissected_mouse/03_concatenated_data|breast_cancer/data/dissected_mouse/')
tissue_labels = c(colon = 'Colon', bone_marrow = 'Bone marrow', fetal_liver = 'Fetal liver', breast_cancer = 'Breast cancer')
score_suffix = '_MAP_SAS'      # benchmark definition whose unintegrated scores are shown ('' = reference definition)
species_cols = c(human = '#4477AA', mouse = '#EE6677')    # as analysis/visualisations.R
out = 'figures/benchmark_summary/unintegrated_umap'
tmp_dir = 'figures/benchmark_summary/unintegrated_umap_data'
n_per_species = Inf            # cells drawn per species and tissue (all if fewer); Inf = every evaluation cell
zoom = c(0.01, 0.99)           # axis limits: central 98% of the drawn cells
point_size = 0.55
point_alpha = 1
# distinct colours for up to 15 cell types (Tol bright / muted, Okabe-Ito); assigned in sorted label order
ct_palette = c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#E69F00', '#332288', '#8C564B',
               '#BBBBBB', '#44AA99', '#000000', '#9467BD', '#D62728', '#882255')   # no two yellows / pinks in one tissue

dir.create(tmp_dir, showWarnings = FALSE)
read_case = function(tissue) {
  parts = strsplit(cases[[tissue]], '|', fixed = TRUE)[[1]]
  csv = file.path(tmp_dir, paste0(tissue, '.csv'))
  if (!file.exists(csv)) stop('run analysis/appendix_unintegrated_umap_compute.py first: ', csv)
  d = read_csv(csv, show_col_types = FALSE)
  eval_cells = readLines(file.path(parts[2], 'eval_cells.txt'))
  d = d[d$cell %in% eval_cells, ]
  stopifnot(nrow(d) == length(eval_cells))
  d
}
sc = read_csv(sprintf('figures/benchmark_summary/benchmark_improvement%s.csv', score_suffix), show_col_types = FALSE) |>
  distinct(tissue, bio = ref_bio, batch = ref_batch, overall = ref_overall)
stopifnot(nrow(sc) == length(cases))
write_csv(sc, paste0(out, '_scores.csv'))

set.seed(1)
draw_check = tibble()
counts = tibble()
row_plot = function(tissue) {
  full = read_case(tissue) |> mutate(species = sub('^hsapiens$', 'human', sub('^mmusculus$', 'mouse', species)))
  d = (if (is.finite(n_per_species)) full |> group_by(species) |> slice_sample(n = n_per_species) |> ungroup() else full) |>
    slice_sample(prop = 1)   # random draw order: no species on top
  # check: the cells drawn last (on top) have the same species mix as all drawn cells
  last = tail(d$species, ceiling(nrow(d) / 10))
  draw_check <<- bind_rows(draw_check, tibble(tissue = tissue, mouse_all = mean(d$species == 'mouse'), mouse_last10pct = mean(last == 'mouse')))
  s = sc[sc$tissue == tissue, ]
  xl = quantile(d$umap_1, zoom); yl = quantile(d$umap_2, zoom)
  xl = xl + c(-0.04, 0.04) * diff(xl); yl = yl + c(-0.04, 0.04) * diff(yl)
  types = sort(unique(d$celltype)); stopifnot(length(types) <= length(ct_palette))
  ct_cols = setNames(ct_palette[seq_along(types)], types)
  panel = function(col, values, name) ggplot(d, aes(umap_1, umap_2, colour = .data[[col]])) +
    ggrastr::rasterise(geom_point(size = point_size, stroke = 0, alpha = point_alpha), dpi = 300) +   # points as an image in the PDF, text stays vector
    scale_colour_manual(values = values, name = name) +
    guides(colour = guide_legend(override.aes = list(size = 3.5, alpha = 1), ncol = if (col == 'celltype' && length(values) > 12) 2 else 1)) +
    coord_equal(xlim = xl, ylim = yl) +
    theme_void(base_size = 12) + theme(panel.border = element_rect(fill = NA, colour = 'grey80'),
                                       legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = 'bold'))
  n_sp = table(full$species)
  counts <<- bind_rows(counts, tibble(tissue = tissue, human = n_sp[['human']], mouse = n_sp[['mouse']]))
  # tissue name and its unintegrated scores in their own strip above the panels, so they never overlap the points
  header = ggplot() + theme_void() +
    labs(title = tissue_labels[[tissue]], subtitle = sprintf('bio %.3f · batch %.3f · overall %.3f', s$bio, s$batch, s$overall)) +
    theme(plot.title = element_text(face = 'bold', size = 14), plot.subtitle = element_text(size = 12, colour = 'grey25'))
  wrap_elements(full = header / (panel('species', species_cols, 'Species') | panel('celltype', ct_cols, 'Cell type')) +
                  plot_layout(heights = c(0.001, 1)))
}
p = wrap_plots(lapply(names(cases), row_plot), ncol = 1)
ggsave(paste0(out, '.png'), p, width = 13, height = 22, dpi = 300, bg = 'white')
ggsave(paste0(out, '.pdf'), p, width = 13, height = 22, device = cairo_pdf, bg = 'white')
writeLines(paste('UMAP of the unintegrated data (log-normalised joint PCA of human and mouse, no integration; the reference of the benchmark),',
                 paste('dissected mouse data, evaluation cells; left: species, right: cell type.',
                       if (is.finite(n_per_species)) sprintf(paste('For visibility, %s cells per species are drawn (random sample, so both',
                                                                  'species are equally represented).'), format(n_per_species, big.mark = ','))
                       else 'All evaluation cells are drawn, in random order.',
                       'Axes zoomed to the central 98% of cells.',
                       paste('UMAP settings (as the bone-marrow method figure):', readLines(file.path(tmp_dir, 'settings.txt'))[1])),
                 sprintf(paste('Scores of the unintegrated data in each tissue: bio-conservation = mean of NMI, ARI, ASW label, cLISI%s;',
                               'batch correction = mean of PCR batch, ASW batch, iLISI, graph connectivity%s; overall = 0.6 bio + 0.4 batch.'),
                         if (score_suffix == '_MAP_SAS') ', MAP' else '', if (score_suffix == '_MAP_SAS') ', SAS' else ''),
                 paste(with(inner_join(counts, sc, by = 'tissue'),
                            sprintf('%s: bio %.3f, batch %.3f, overall %.3f; %s human and %s mouse cells.', tissue_labels[tissue], bio, batch, overall,
                                    prettyNum(human, big.mark = ','), prettyNum(mouse, big.mark = ','))), collapse = ' '), sep = '\n'),
           paste0(out, '_caption.txt'))
print(as.data.frame(sc)); print(as.data.frame(draw_check))
