# Marker AUROC of every annotator, all tissues and species in one figure (tissue rows x species columns).
# Values: figures/<tissue>/annotation/marker_auroc.csv (analysis/visualisations.R: AUROC of the mean panel-marker expression for
# the cells given a label vs all other cells; the per-dataset heatmaps marker_auroc_<species>.png show the same numbers).
# As those heatmaps: cell types no annotator calls are dropped; Mean = each annotator's mean over the types it calls.
# style = 'heatmap': tiles as the per-dataset heatmaps (plot_marker_auroc in analysis/visualisations.R: annotators in columns, value
#   printed, diverging fill around 0.5, '-' = no markers or never called). style = 'dots': x = AUROC, annotators offset
#   within each row; colour = method, filled = same-species reference, open = other species' reference.
#   Rscript --vanilla analysis/marker_auroc_figure.R
# Output: figures/cross_tissue/marker_auroc_all.{png,pdf,csv}
suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2); library(patchwork)})

# ---------------- config ----------------
tissue_labels  = c(colon = 'Colon', bone_marrow = 'Bone marrow', fetal_liver = 'Fetal liver', breast_cancer = 'Breast cancer')
species_labels = c(human = 'Human', mouse = 'Mouse')
annotators = c(seurat = 'Seurat', seurat_cross = 'Seurat (cross-species)', singleR = 'SingleR',
               singleR_cross = 'SingleR (cross-species)', ebcc = 'EBCC consensus')
ann_colour = c(seurat = '#228833', seurat_cross = '#228833', singleR = '#AA3377', singleR_cross = '#AA3377', ebcc = '#222222')
ann_fill   = c(seurat = '#228833', seurat_cross = 'white',   singleR = '#AA3377', singleR_cross = 'white',   ebcc = '#222222')
ann_shape  = c(seurat = 21, seurat_cross = 21, singleR = 21, singleR_cross = 21, ebcc = 23)
offset     = c(seurat = 0.28, seurat_cross = 0.14, singleR = 0, singleR_cross = -0.14, ebcc = -0.28)   # within a row
mean_row   = 'Mean (all called types)'
x_lim      = c(0, 1)
style      = 'heatmap'      # 'heatmap' (as marker_auroc_<species>.png) or 'dots'
out        = paste0('figures/cross_tissue/marker_auroc_all', if (style == 'dots') '_dots' else '')
# -----------------------------------------

d = bind_rows(lapply(names(tissue_labels), \(t)
  read_csv(sprintf('figures/%s/annotation/marker_auroc.csv', t), show_col_types = FALSE) |> mutate(tissue = t))) |>
  group_by(tissue, species, celltype) |> filter(any(!is.na(auroc))) |> ungroup()    # as plot_marker_auroc
stopifnot(setequal(unique(d$annotator), names(annotators)), all(d$auroc >= 0 & d$auroc <= 1, na.rm = TRUE))
means = d |> group_by(tissue, species, annotator) |> summarise(auroc = mean(auroc, na.rm = TRUE), .groups = 'drop') |>
  mutate(celltype = mean_row)
dd = bind_rows(d, means) |> mutate(row_lab = ifelse(celltype == mean_row, mean_row, sprintf('%s (%d)', gsub('_', ' ', celltype), n_markers)))   # mouse names: '_' shown as spaces
write_csv(dd |> select(tissue, species, annotator, celltype, n_markers, n_cells, auroc), paste0(out, '.csv'))

panel = function(t, s, show_legend = FALSE) {
  x = filter(dd, tissue == t, species == s)
  rows = c(unique(x$row_lab[x$celltype != mean_row]), mean_row)          # marker order, mean at the bottom
  x = x |> mutate(y = match(row_lab, rev(rows)) + offset[annotator], annotator = factor(annotator, levels = names(annotators)))
  ggplot(x, aes(auroc, y)) +
    annotate('rect', xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.5, fill = 'grey94') +             # mean row band
    geom_hline(yintercept = seq_along(rows), colour = 'grey90', linewidth = 0.3) +
    geom_vline(xintercept = 0.5, linetype = 'dashed', colour = 'grey45', linewidth = 0.4) +
    geom_point(aes(colour = annotator, fill = annotator, shape = annotator), size = 2.1, stroke = 0.7, na.rm = TRUE) +
    scale_colour_manual(values = ann_colour, labels = annotators, name = NULL, drop = FALSE) +
    scale_fill_manual(values = ann_fill, labels = annotators, name = NULL, drop = FALSE) +
    scale_shape_manual(values = ann_shape, labels = annotators, name = NULL, drop = FALSE) +
    scale_y_continuous(breaks = seq_along(rows), labels = rev(rows), expand = expansion(add = 0.4)) +
    scale_x_continuous(limits = x_lim, breaks = seq(0, 1, 0.25), labels = c('0', '0.25', '0.5', '0.75', '1')) +
    labs(x = 'Marker AUROC', y = NULL, title = sprintf('%s, %s', tissue_labels[[t]], tolower(species_labels[[s]]))) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = 'bold', size = 10),
          axis.text.y = element_text(size = 8, face = ifelse(rev(rows) == mean_row, 'bold', 'plain')),
          legend.position = if (show_legend) 'bottom' else 'none')
}
heat_panel = function(t, s) {   # as plot_marker_auroc (analysis/visualisations.R)
  x = filter(dd, tissue == t, species == s)
  rows = c(unique(x$row_lab[x$celltype != mean_row]), mean_row)
  x = x |> mutate(row = factor(row_lab, levels = rev(rows)), annotator = factor(annotator, levels = heat_cols))
  ggplot(x, aes(annotator, row, fill = auroc)) +
    geom_tile(colour = 'white') +
    geom_hline(yintercept = 1.5, linewidth = 0.6) +
    geom_text(aes(label = ifelse(is.na(auroc), '\u2013', sprintf('%.2f', auroc)),
                  fontface = ifelse(celltype == mean_row, 'bold', 'plain')), size = 2.9) +
    scale_fill_gradient2(low = '#B2182B', mid = 'white', high = '#2166AC', midpoint = 0.5, limits = c(0.3, 1),
                         oob = scales::squish, na.value = 'grey90', name = 'Marker AUROC') +
    scale_x_discrete(labels = heat_lab, drop = FALSE) +
    labs(x = NULL, y = NULL, title = sprintf('%s, %s', tissue_labels[[t]], tolower(species_labels[[s]]))) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(),
          plot.title = element_text(face = 'bold', size = 10))
}
heat_cols = c('seurat', 'singleR', 'seurat_cross', 'singleR_cross', 'ebcc')   # column order of the per-dataset heatmaps
heat_lab  = c(seurat = 'Seurat', singleR = 'SingleR', seurat_cross = 'Seurat\n(cross)', singleR_cross = 'SingleR\n(cross)', ebcc = 'EBCC')
if (style == 'heatmap') panel = function(t, s, show_legend = TRUE) heat_panel(t, s)
n_rows = sapply(names(tissue_labels), \(t) max(sapply(names(species_labels), \(s)
  n_distinct(filter(dd, tissue == t, species == s)$row_lab))))
ps = unlist(lapply(names(tissue_labels), \(t) lapply(names(species_labels), \(s) panel(t, s, show_legend = TRUE))), recursive = FALSE)
p = wrap_plots(ps, ncol = 2) + plot_layout(heights = n_rows + 2.5, guides = 'collect') + plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'bottom', plot.tag = element_text(face = 'bold', size = 12))
h = 1.2 + (if (style == 'heatmap') 0.3 else 0.32) * sum(n_rows + 2.5)
w = if (style == 'heatmap') 11 else 12
ggsave(paste0(out, '.png'), p, width = w, height = h, dpi = 300, bg = 'white')
ggsave(paste0(out, '.pdf'), p, width = w, height = h, bg = 'white', device = cairo_pdf)
print(as.data.frame(means |> mutate(auroc = round(auroc, 3)) |> tidyr::pivot_wider(names_from = annotator, values_from = auroc) |> select(-celltype)), row.names = FALSE)
cat('wrote', paste0(out, '.{png,pdf,csv}'), sprintf('(%.1f in tall)\n', h))
