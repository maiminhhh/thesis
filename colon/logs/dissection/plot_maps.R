setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
suppressMessages(source('../library.R'))
d = readRDS('logs/dissection/mouse_cells.rds'); cells = d$cells
print(summary(cells[c('x', 'y')]))
print(sapply(cells[names(d$offtarget)], quantile, c(.5, .9, .99, .999)), digits = 2)
bins = cells |> group_by(bx, by) |> summarise(n = n(),
  frac_epithelial = mean(grepl('^Epithelial', label)), frac_smooth_muscle = mean(label == 'Smooth muscle'),
  frac_immune = mean(label %in% c('Immune', 'B/Plasma')),
  across(all_of(names(d$offtarget)), mean), .groups = 'drop') |> filter(n >= 20)
saveRDS(bins, 'logs/dissection/mouse_bins.rds')
long = bins |> pivot_longer(c(starts_with('frac_'), all_of(names(d$offtarget))), names_to = 'metric')
# ImageDimPlot orientation: y horizontal, x vertical
p = ggplot(long, aes(by + d$B / 2, bx + d$B / 2, fill = value)) + geom_tile() +
  facet_wrap(~ metric, ncol = 3) + scale_fill_viridis_c() + coord_fixed() + scale_y_reverse() +
  annotate('rect', xmin = 2590, xmax = 4625, ymin = 260, ymax = 2685, fill = NA, colour = 'red') +
  labs(x = 'y', y = 'x', title = 'mouse colon section, 250 x 250 bins (red: current cut)') + theme_bw()
ggsave('logs/dissection/bin_maps.png', p, width = 14, height = 12, dpi = 110)
p2 = ggplot(cells, aes(y, x, colour = label)) + geom_point(size = 0.05) + scale_y_reverse() + coord_fixed() +
  annotate('rect', xmin = 2590, xmax = 4625, ymin = 260, ymax = 2685, fill = NA, colour = 'black') +
  guides(colour = guide_legend(override.aes = list(size = 3))) + theme_bw() + labs(title = 'EBCC labels')
ggsave('logs/dissection/labels.png', p2, width = 12, height = 9, dpi = 110)
