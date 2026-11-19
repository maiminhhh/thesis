# Dissection record (../EBCC_PIPELINE.md section 5b, step 4) with the final mouse EBCC labels:
# composition before / after the cut plus the human composition, and the cut over the whole section.
# The cut itself is applied in 03_concatenate.qmd (pipeline/03_partA_datasets.R); keep the two in sync.
suppressMessages({source('../library.R'); source('../source.R')})
cut = list(name = 'central tumour block', xmin = 4500, xmax = 8500, ymin = 1500, ymax = 5500)

mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
human = readRDS('data/02_annotated_data/02_human_annotated.rds')
coords = GetTissueCoordinates(mouse)
cells = tibble(cell = coords$cell, x = coords$x, y = coords$y,
               label = as.character(mouse@meta.data[coords$cell, 'ebcc_label'])) |>
  mutate(in_cut = x >= cut$xmin & x <= cut$xmax & y >= cut$ymin & y <= cut$ymax)

pct = function(l) { t = table(l); setNames(as.numeric(100 * t / sum(t)), names(t)) }
labs_all = sort(unique(c(cells$label, as.character(human$ebcc_label))))
compo = tibble(label = labs_all,
               mouse_whole_section_pct = pct(cells$label)[labs_all],
               mouse_dissected_pct     = pct(cells$label[cells$in_cut])[labs_all],
               human_pct               = pct(as.character(human$ebcc_label))[labs_all]) |>
  mutate(across(ends_with('_pct'), \(v) replace_na(v, 0)))
compo = bind_rows(compo, tibble(label = 'n_cells', mouse_whole_section_pct = nrow(cells),
                                mouse_dissected_pct = sum(cells$in_cut), human_pct = ncol(human)))
write_csv(compo, 'data/02_annotated_data/mouse_dissection_composition.csv')
print(as.data.frame(compo |> mutate(across(where(is.numeric), \(v) round(v, 1)))), row.names = FALSE)

# plotted as in ImageDimPlot: y horizontal, x vertical
p = ggplot(cells, aes(y, x, colour = label)) +
  geom_point(size = 0.01) +
  annotate('rect', xmin = cut$ymin, xmax = cut$ymax, ymin = cut$xmin, ymax = cut$xmax,
           fill = NA, colour = 'black', linewidth = 0.8) +
  coord_fixed() + guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = 'y (horizontal in ImageDimPlot)', y = 'x (vertical in ImageDimPlot)', colour = 'EBCC label',
       title = sprintf('breast_cancer mouse dissection: %s', cut$name),
       subtitle = sprintf('%d <= x <= %d, %d <= y <= %d; %d / %d cells', cut$xmin, cut$xmax, cut$ymin, cut$ymax,
                          sum(cells$in_cut), nrow(cells))) +
  theme_minimal()
dir.create('../figures/breast_cancer/spatial', showWarnings = FALSE, recursive = TRUE)
ggsave('../figures/breast_cancer/spatial/mouse_dissection.png', p, width = 9, height = 9, bg = 'white')
cat('dissection record done\n')
