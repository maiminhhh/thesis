# Combined marker dot plot: one panel per tissue (rows) and species (columns), cells grouped by EBCC label.
# Values from analysis/marker_dotplot_compute.py (as Seurat's DotPlot): dot size = % of cells with the gene detected; colour =
# mean expression (mean of expm1 of the log-normalised value, log1p), scaled across the labels of the panel per gene and
# clipped at +-2.5, as DotPlot's avg.exp.scaled. Genes grouped by the cell type they mark (strips), at most
# max_per_type per cell type (the first listed in the annotation script; full lists in marker_table.tex).
#   Rscript --vanilla analysis/marker_dotplot_figure.R
#   ... analysis/marker_dotplot_figure.R selection   -> the EBCC frac-selection half only (analysis/marker_dotplot_compute.py selection)
# Output: figures/cross_tissue/marker_dotplot_ebcc[_<half>].{png,pdf}
suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2); library(patchwork)})

# ---------------- config ----------------
half         = if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else 'all'   # all | selection | evaluation
sfx          = if (half == 'all') '' else paste0('_', half)
in_file      = sprintf('figures/cross_tissue/marker_dotplot_values%s.csv', sfx)
markers_file = 'figures/cross_tissue/marker_table.csv'
out          = paste0('figures/cross_tissue/marker_dotplot_ebcc', sfx)
max_per_type = 5            # markers shown per cell type (Inf = all)
clip         = 2.5          # as DotPlot (col.min / col.max)
tissue_labels = c(colon = 'Colon', bone_marrow = 'Bone marrow', fetal_liver = 'Fetal liver', breast_cancer = 'Breast cancer')
species_labels = c(human = 'Human', mouse = 'Mouse')
col_low = '#eef1f6'; col_high = '#08306b'   # one-hue sequential (light = low)
# -----------------------------------------

v = read_csv(in_file, show_col_types = FALSE)
mt = read_csv(markers_file, show_col_types = FALSE)

panel = function(t, s) {
  x = filter(v, tissue == t, species == s)
  grp_order = filter(mt, tissue == t, species == s, !is.na(on_panel)) |> pull(label)   # marker groups, table order
  genes = x |> distinct(gene, marker_of) |> mutate(marker_of = factor(marker_of, levels = grp_order)) |>
    left_join(x |> distinct(gene) |> mutate(pos = row_number()), by = 'gene') |>       # order as in the marker list
    arrange(marker_of, pos) |> group_by(marker_of) |> slice_head(n = max_per_type) |> ungroup()
  rows = c(intersect(grp_order, unique(x$label)), sort(setdiff(unique(x$label), grp_order)))
  d = x |> filter(gene %in% genes$gene) |> group_by(gene) |>
    mutate(z = as.numeric(scale(log1p(avg_exp))), z = pmin(pmax(ifelse(is.na(z), 0, z), -clip), clip)) |> ungroup() |>
    mutate(gene = factor(gene, levels = genes$gene), label = factor(label, levels = rev(rows)),
           marker_of = factor(marker_of, levels = grp_order))
  ggplot(d, aes(gene, label)) +
    geom_point(aes(size = pct, colour = z)) +
    facet_grid(. ~ marker_of, scales = 'free_x', space = 'free_x') +
    scale_colour_gradient(low = col_low, high = col_high, limits = c(-clip, clip), name = 'Scaled mean\nexpression') +
    scale_size_area(max_size = 3.2, limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), name = '% of cells\ndetected') +
    labs(x = NULL, y = NULL, title = sprintf('%s, %s%s', tissue_labels[[t]], tolower(species_labels[[s]]),
                                             if (half == 'all') '' else sprintf(' (%s half: %s cells)', half,
                                               scales::comma(sum(distinct(x, label, n_cells)$n_cells))))) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, face = 'italic'),
          axis.text.y = element_text(size = 7), panel.grid.major = element_line(colour = 'grey93', linewidth = 0.3),
          panel.grid.minor = element_blank(), panel.spacing.x = unit(1.5, 'pt'),
          strip.text.x = element_text(angle = 90, hjust = 0, size = 6.5), strip.background = element_rect(fill = 'grey95', colour = 'grey80'),
          plot.title = element_text(size = 9, face = 'bold'))
}
ps = unlist(lapply(names(tissue_labels), \(t) lapply(names(species_labels), \(s) panel(t, s))), recursive = FALSE)
p = wrap_plots(ps, ncol = 2) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'bottom', legend.box = 'horizontal', plot.tag = element_text(face = 'bold', size = 11))
ggsave(paste0(out, '.png'), p, width = 16, height = 20, dpi = 300, bg = 'white')
ggsave(paste0(out, '.pdf'), p, width = 16, height = 20, bg = 'white', device = cairo_pdf)
n = v |> distinct(tissue, species, gene, marker_of) |> group_by(tissue, species, marker_of) |>
  summarise(k = n(), .groups = 'drop') |> mutate(shown = pmin(k, max_per_type)) |> group_by(tissue, species) |>
  summarise(markers = sum(k), shown = sum(shown), .groups = 'drop')
print(as.data.frame(n), row.names = FALSE)
cat('wrote', paste0(out, '.{png,pdf}'), '\n')
