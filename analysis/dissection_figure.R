# Figure: how the mouse sections were dissected.  One row per tissue:
#   A  the full mouse section coloured by EBCC label, the dissected region outlined (cells outside it faded);
#   B  the human sample, coloured the same way;
#   (the composition of the human sample, the dissected region and the full section, and the off-target labels named
#   in EBCC_PIPELINE.md section 5b, are written to mouse_dissection_composition.csv and printed, not plotted).
# Input: figures/dissection/dissection_cells.csv.gz (analysis/dissection_figure_compute.py: raw centroids + EBCC labels; the
# dissection recomputed from each tissue's rectangle reproduces the saved dissected cells).
# Plots x horizontally and y vertically (y increasing downwards, as in the images). Labels below 1% in all three groups
# of a tissue are pooled as 'Other'.
#   Rscript --vanilla analysis/dissection_figure.R
# Output: figures/dissection/mouse_dissection.{png,pdf} (point layers rasterised), mouse_dissection_composition.csv
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)})

# ---------------- config ----------------
in_file  = 'figures/dissection/dissection_cells.csv.gz'
out      = 'figures/dissection'
tissues  = c(colon = 'Colon', bone_marrow = 'Bone marrow', fetal_liver = 'Fetal liver', breast_cancer = 'Breast cancer')
rects    = list(colon = c(260, 2685, 2590, 4625), bone_marrow = c(500, 4750, 2500, 3500),   # xmin, xmax, ymin, ymax
                fetal_liver = c(0, 3100, 500, 3000), breast_cancer = c(4500, 8500, 1500, 5500))
off_target = list(bone_marrow = c('MSC', 'T cell', 'NK/ILC'),           # cartilage / cortical stroma (EBCC_PIPELINE.md 5b)
                  breast_cancer = c('Normal Epithelial'))               # normal mammary gland left out (03_concatenate.qmd)
min_pct  = 1                                                           # labels below this in all groups -> 'Other'
cols     = c('#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499',
             '#EE8866', '#6699CC', '#BBCC33')                          # Paul Tol muted + light; fixed order by human abundance
other_col = 'grey80'
pt_fill  = 1.1                                                         # point diameter / typical spacing between cells
panel_in = 4.6                                                         # approximate panel width (inches) for point sizing
raster_dev = 'ragg_png'                                                # ragg keeps ggplot's point sizes (the default cairo device enlarges them)
raster_dpi = 400                                                       # point layers are rasterised (PDF stays small); the rest is vector
fade     = 0.35                                                        # alpha of mouse cells outside the dissected region
# -----------------------------------------

d = read_csv(in_file, show_col_types = FALSE, col_types = cols(dissected = 'c')) |>
  mutate(dissected = dissected %in% c('True', 'TRUE', '1.0', '1'))
groups = c(human = 'Human sample', dissected = 'Mouse, dissected', full = 'Mouse, full section')
comp = bind_rows(d |> filter(species == 'human') |> mutate(group = 'human'),
                 d |> filter(species == 'mouse', dissected) |> mutate(group = 'dissected'),
                 d |> filter(species == 'mouse') |> mutate(group = 'full')) |>
  count(tissue, group, label) |> group_by(tissue, group) |> mutate(pct = 100 * n / sum(n)) |> ungroup()
keep = comp |> group_by(tissue, label) |> summarise(keep = max(pct) >= min_pct, .groups = 'drop')
comp = comp |> left_join(keep, by = c('tissue', 'label')) |> mutate(label2 = ifelse(keep, label, 'Other')) |>
  group_by(tissue, group, label2) |> summarise(n = sum(n), pct = sum(pct), .groups = 'drop')
write_csv(comp |> rename(label = label2) |> mutate(pct = round(pct, 2)) |> arrange(tissue, group, -pct),
          file.path(out, 'mouse_dissection_composition.csv'))
d = d |> left_join(keep, by = c('tissue', 'label')) |> mutate(label2 = ifelse(keep, label, 'Other'))

# total variation distance to the human composition (printout and caption)
tv = comp |> select(tissue, group, label2, pct) |> pivot_wider(names_from = group, values_from = pct, values_fill = 0) |>
  group_by(tissue) |> summarise(full = sum(abs(full - human)) / 200, dissected = sum(abs(dissected - human)) / 200)
cat('Total variation distance to the human composition (0 = identical):\n'); print(as.data.frame(tv), digits = 3)
cat('Off-target labels (% of cells, full section -> dissected):\n')
for (t in names(off_target)) comp |> filter(tissue == t, label2 %in% off_target[[t]], group != 'human') |>
  select(label2, group, pct) |> pivot_wider(names_from = group, values_from = pct) |>
  with(cat(sprintf('  %-14s %-18s %5.1f -> %5.1f\n', t, label2, full, dissected), sep = ''))

theme_map = theme_void(base_size = 10) +
  theme(plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 3)), legend.position = 'none')
row_plot = function(t) {
  lev = comp |> filter(tissue == t, group == 'human', label2 != 'Other') |> arrange(-pct) |> pull(label2)
  lev = c(lev, setdiff(unique(comp$label2[comp$tissue == t & comp$label2 != 'Other']), lev))   # mouse-only labels last
  pal = setNames(cols[seq_along(lev)], lev)
  if (any(comp$label2[comp$tissue == t] == 'Other')) { lev = c(lev, 'Other'); pal = c(pal, Other = other_col) }
  x = d |> filter(tissue == t) |> mutate(label2 = factor(label2, levels = lev))
  m = filter(x, species == 'mouse'); h = filter(x, species == 'human'); r = rects[[t]]
  fmt = \(v) format(v, big.mark = ',')
  # point size from cell density: typical spacing sqrt(area / n), in data units, scaled to the panel width (mm / .pt)
  psize = function(z) { w = diff(range(z$x)); a = diff(range(z$x)) * diff(range(z$y))
    min(0.6, max(0.03, pt_fill * sqrt(a / nrow(z)) / w * panel_in * 25.4 / .pt * 2)) }
  sm = psize(m); sh = psize(h)
  pA = ggplot(m, aes(x, y, colour = label2)) +
    ggrastr::rasterise(geom_point(data = filter(m, !dissected), size = sm, stroke = 0, alpha = fade), dpi = raster_dpi, dev = raster_dev) +
    ggrastr::rasterise(geom_point(data = filter(m, dissected), size = sm, stroke = 0), dpi = raster_dpi, dev = raster_dev) +
    annotate('rect', xmin = r[1], xmax = r[2], ymin = r[3], ymax = r[4], fill = NA, colour = 'black', linewidth = 0.6) +
    scale_colour_manual(values = pal, drop = FALSE) + scale_y_reverse() + coord_fixed() +
    labs(title = sprintf('Mouse section (%s cells; outlined: dissected, %s)', fmt(nrow(m)), fmt(sum(m$dissected)))) + theme_map
  pB = ggplot(h, aes(x, y, colour = label2)) +
    ggrastr::rasterise(geom_point(size = sh, stroke = 0), dpi = raster_dpi, dev = raster_dev) +
    scale_colour_manual(values = pal, drop = FALSE) + scale_y_reverse() + coord_fixed() +
    labs(title = sprintf('Human sample (%s cells)', fmt(nrow(h)))) + theme_map +
    guides(colour = guide_legend(override.aes = list(size = 2.5), ncol = 1)) +
    theme(legend.position = 'right', legend.title = element_blank(), legend.text = element_text(size = 8.5),
          legend.key.size = unit(0.35, 'cm'))
  lab = wrap_elements(full = grid::textGrob(tissues[[t]], rot = 90, gp = grid::gpar(fontsize = 12, fontface = 'bold')))
  cat(sprintf('%-14s point size mouse %.3f, human %.3f\n', t, sm, sh))
  (lab | pA | pB) + plot_layout(widths = c(0.08, 1.3, 1))
}
p = wrap_plots(lapply(names(tissues), row_plot), ncol = 1) + plot_annotation(tag_levels = list(c('', 'A', 'B', rep('', 12))))
ggsave(file.path(out, 'mouse_dissection.png'), p, width = 12, height = 15, dpi = 300, bg = 'white')
ggsave(file.path(out, 'mouse_dissection.pdf'), p, width = 12, height = 15, bg = 'white', device = cairo_pdf)
cat('wrote', file.path(out, 'mouse_dissection.{png,pdf}'), '\n')
