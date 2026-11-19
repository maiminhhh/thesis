# Mouse dissection map (../EBCC_PIPELINE.md section 5b) for breast_cancer.
# Human sample: a breast tumour TMA core (MERFISH, HubcTMA region R4), so the matching mouse compartment is the
# tumour mass (tumour epithelium + intratumoural stroma), without normal mammary gland, fat, lymph node or skin.
# Bins GetTissueCoordinates() into 250 x 250 squares; per bin: fraction of cells whose label is Cancer Epithelial
# (target) and mean expression of off-target marker sets (label-free). Optional candidate rectangles are drawn
# and their composition is compared with the whole section and the human sample.
#   Rscript pipeline/mouse_dissection_map.R <label column> <output tag> [candidates.csv]
# candidates.csv: name,xmin,xmax,ymin,ymax (GetTissueCoordinates units; in ImageDimPlot y is horizontal, x vertical)
suppressMessages({source('../library.R'); source('../source.R')})
args = commandArgs(trailingOnly = TRUE)
label_col = args[1]; tag = args[2]; cand_file = if (length(args) > 2) args[3] else NA
bin = 250
fig_dir = '../figures/breast_cancer/spatial'
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
coords = GetTissueCoordinates(mouse)                      # columns x, y, cell
lab = as.character(mouse@meta.data[coords$cell, label_col])
X = LayerData(mouse, assay = DefaultAssay(mouse), layer = 'data')
sets = list(tumour = c('Mki67', 'Top2a', 'Cdh1', 'Cldn4'),               # reference only (on-target)
            normal_basal = c('Krt5', 'Krt14', 'Krt17', 'Trp63'),         # normal mammary gland
            adipose = c('Lep', 'Pparg'),                                 # fat
            lymph_node = c('Ccr7', 'Cd3g'))                              # lymph node
sets = lapply(sets, intersect, rownames(X))
score = sapply(sets, function(g) Matrix::colMeans(X[g, coords$cell, drop = FALSE]))

cells = tibble(cell = coords$cell, x = coords$x, y = coords$y, label = lab) |>
  bind_cols(as_tibble(score)) |>
  mutate(bx = floor(x / bin), by = floor(y / bin))
bins = cells |>
  dplyr::group_by(bx, by) |>
  summarise(n = n(), frac_cancer_epithelial = mean(label %in% 'Cancer Epithelial'),
            frac_normal_epithelial = mean(label %in% 'Normal Epithelial'),
            across(all_of(names(sets)), mean), .groups = 'drop') |>
  filter(n >= 20)
write_csv(bins, sprintf('data/02_annotated_data/mouse_dissection_bins_%s.csv', tag))

cands = if (!is.na(cand_file)) read_csv(cand_file, show_col_types = FALSE) else
  tibble(name = 'current', xmin = 8560, xmax = 11745, ymin = 920, ymax = 5455)
# plotted as in ImageDimPlot: y horizontal, x vertical
rects = geom_rect(data = cands, aes(xmin = ymin, xmax = ymax, ymin = xmin, ymax = xmax),
                  inherit.aes = FALSE, fill = NA, colour = 'red', linewidth = 0.6)
rect_labels = geom_text(data = cands, aes(x = ymin, y = xmax, label = name), inherit.aes = FALSE,
                        hjust = 0, vjust = -0.3, colour = 'red', size = 3)
long = bins |>
  pivot_longer(c(frac_cancer_epithelial, frac_normal_epithelial, all_of(names(sets))), names_to = 'layer') |>
  group_by(layer) |> mutate(value_scaled = value / max(value, na.rm = TRUE)) |> ungroup()
p_map = ggplot(long, aes(by * bin + bin / 2, bx * bin + bin / 2, fill = value_scaled)) +
  geom_tile(width = bin, height = bin) + rects + rect_labels +
  facet_wrap(~ layer, nrow = 2) + coord_fixed() + scale_fill_viridis_c(name = 'scaled\n(max = 1)') +
  labs(x = 'y (horizontal in ImageDimPlot)', y = 'x (vertical in ImageDimPlot)',
       title = sprintf('breast_cancer mouse section, %d x %d bins, labels: %s', bin, bin, label_col),
       subtitle = 'fractions use the labels; normal_basal / adipose / lymph_node / tumour are mean marker expression') +
  theme_minimal()
ggsave(file.path(fig_dir, sprintf('mouse_dissection_map_%s.png', tag)), p_map, width = 14, height = 9, bg = 'white')
p_lab = ggplot(cells, aes(y, x, colour = label)) + geom_point(size = 0.01) + rects + rect_labels +
  coord_fixed() + guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = 'y (horizontal in ImageDimPlot)', y = 'x (vertical in ImageDimPlot)', colour = label_col,
       title = 'mouse section by label, with candidate cuts') + theme_minimal()
ggsave(file.path(fig_dir, sprintf('mouse_dissection_labels_%s.png', tag)), p_lab, width = 9, height = 9, bg = 'white')

# composition (label %) of the whole section, each candidate, and the human sample (EBCC)
comp = function(l) tibble(label = l) |> dplyr::count(label) |> mutate(pct = 100 * n / sum(n))
inside = function(r) cells$x >= r$xmin & cells$x <= r$xmax & cells$y >= r$ymin & cells$y <= r$ymax
human = readRDS('data/02_annotated_data/02_human_annotated.rds')
compo = bind_rows(
  comp(cells$label) |> mutate(set = 'mouse: whole section'),
  map_dfr(seq_len(nrow(cands)), function(i) comp(cells$label[inside(cands[i, ])]) |>
            mutate(set = sprintf('mouse: %s', cands$name[i]))),
  comp(as.character(human$ebcc_label)) |> mutate(set = 'human (EBCC)'))
write_csv(compo, sprintf('data/02_annotated_data/mouse_dissection_options_%s.csv', tag))
wide = compo |> select(set, label, pct) |> pivot_wider(names_from = set, values_from = pct, values_fill = 0)
print(as.data.frame(wide |> mutate(across(where(is.numeric), \(v) round(v, 1)))))
cat('cells:', paste(sprintf('%s = %d', c('whole', cands$name), c(nrow(cells), sapply(seq_len(nrow(cands)), \(i) sum(inside(cands[i, ]))))), collapse = ', '), '\n')
