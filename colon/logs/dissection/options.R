setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
suppressMessages(source('../library.R'))
d = readRDS('logs/dissection/mouse_cells.rds'); cells = d$cells; B = d$B
opts = tribble(~option, ~xmin, ~xmax, ~ymin, ~ymax,
  'A current cut',          260, 2685, 2590, 4625,
  'B upper colon segment',    0, 2100,    0, 7300,
  'C right colon loop',    3000, 4950, 5100, 7150,
  'D whole section',          0, 6100,    0, 7300)
inside = function(o) cells$x >= o$xmin & cells$x <= o$xmax & cells$y >= o$ymin & cells$y <= o$ymax
comp = bind_rows(lapply(seq_len(nrow(opts)), function(i) {
  o = opts[i, ]; k = inside(o)
  tibble(option = o$option, n = sum(k), label = names(table(cells$label[k])), pct = 100 * as.numeric(prop.table(table(cells$label[k]))))
}))
human = read.csv(text = system(paste0(file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'), ' -c "import anndata as ad,warnings;warnings.filterwarnings(\'ignore\');print(ad.read_h5ad(\'data/02_annotated_data/02_human_annotated.h5ad\',backed=\'r\').obs[\'ebcc_label\'].value_counts().to_csv())"'), intern = TRUE))
names(human) = c('label', 'count'); human$pct = 100 * human$count / sum(human$count)
wide = comp |> select(option, label, pct) |> pivot_wider(names_from = option, values_from = pct, values_fill = 0) |>
  full_join(select(human, label, human = pct), by = 'label') |> arrange(-human)
print(as.data.frame(wide), digits = 2)
print(comp |> distinct(option, n))
off = bind_rows(lapply(seq_len(nrow(opts)), function(i) { k = inside(opts[i, ])
  tibble(option = opts$option[i], adipose_cells_pct = 100 * mean(cells$adipose[k] > 2), squamous_mean = mean(cells$squamous[k]),
         skeletal_mean = mean(cells$skeletal[k]), cartilage_bone_pct = 100 * mean(cells$cartilage_bone[k] > 1), urothelium_pct = 100 * mean(cells$urothelium[k] > 1)) }))
print(as.data.frame(off), digits = 2)
saveRDS(list(opts = opts, comp = comp, human = human, off = off), 'logs/dissection/options.rds')
# maps: each metric on its own colour scale (0 .. its 99th percentile over bins), candidate cuts drawn
bins = readRDS('logs/dissection/mouse_bins.rds')
long = bins |> pivot_longer(c(starts_with('frac_'), all_of(names(d$offtarget))), names_to = 'metric') |>
  group_by(metric) |> mutate(scaled = pmin(1, value / quantile(value, 0.99))) |> ungroup()
p = ggplot(long, aes(by + B / 2, bx + B / 2, fill = scaled)) + geom_tile() +
  geom_rect(data = opts, aes(xmin = ymin, xmax = ymax, ymin = xmin, ymax = xmax, colour = option), inherit.aes = FALSE, fill = NA, linewidth = 0.6) +
  facet_wrap(~ metric, ncol = 3) + scale_fill_viridis_c(name = 'value /\n99th pct') + coord_fixed() + scale_y_reverse() +
  labs(x = 'y (horizontal in ImageDimPlot)', y = 'x (vertical)', title = 'mouse colon: 250 x 250 bins, each metric on its own scale',
       subtitle = 'frac_*: EBCC label fractions; others: mean log-expression of marker sets (label-independent)') + theme_bw()
ggsave('logs/dissection/bin_maps_options.png', p, width = 15, height = 13, dpi = 100)
p2 = ggplot(cells, aes(y, x, colour = label)) + geom_point(size = 0.03) + scale_y_reverse() + coord_fixed() +
  geom_rect(data = opts, aes(xmin = ymin, xmax = ymax, ymin = xmin, ymax = xmax, linetype = option), inherit.aes = FALSE, fill = NA, colour = 'black') +
  guides(colour = guide_legend(override.aes = list(size = 3))) + theme_bw() + labs(title = 'mouse EBCC labels with candidate cuts', x = 'y', y = 'x')
ggsave('logs/dissection/labels_options.png', p2, width = 13, height = 9, dpi = 100)
