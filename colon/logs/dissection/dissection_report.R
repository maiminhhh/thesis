# Section 5b step 4: EBCC label composition of the whole mouse section, the dissection cut (03_concatenate.qmd)
# and the human sample -> data/02_annotated_data/mouse_dissection_composition.csv, and the cut drawn over the
# whole section -> ../figures/colon/spatial/mouse_dissection.png. Uses the final 02 labels.
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
suppressMessages(source('../library.R'))
src = readLines('03_concatenate.qmd')
eval(parse(text = grep('^(xmin|ymin) = ', src, value = TRUE)))   # the cut, as written in 03
mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
coords = GetTissueCoordinates(mouse)
cells = data.frame(x = coords$x, y = coords$y, label = as.character(mouse$ebcc_label[match(coords$cell, colnames(mouse))]))
rm(mouse); gc()
cells$cut = cells$x >= xmin & cells$x <= xmax & cells$y >= ymin & cells$y <= ymax
human = read.csv(text = system(paste0(file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'), ' -c "import anndata as ad,warnings;warnings.filterwarnings(\'ignore\');print(ad.read_h5ad(\'data/02_annotated_data/02_human_annotated.h5ad\',backed=\'r\').obs[\'ebcc_label\'].value_counts().to_csv())"'), intern = TRUE))
names(human) = c('label', 'n')
comp = function(lab, group) { t = table(lab); data.frame(group = group, label = names(t), n = as.integer(t), pct = 100 * as.numeric(t) / sum(t)) }
out = rbind(comp(cells$label, 'mouse_whole_section'), comp(cells$label[cells$cut], 'mouse_dissected'),
            data.frame(group = 'human', label = human$label, n = human$n, pct = 100 * human$n / sum(human$n)))
write.csv(out, 'data/02_annotated_data/mouse_dissection_composition.csv', row.names = FALSE)
print(reshape(out[c('group', 'label', 'pct')], idvar = 'label', timevar = 'group', direction = 'wide'), digits = 2)
p = ggplot(cells, aes(y, x, colour = label)) + geom_point(size = 0.03) + scale_y_reverse() + coord_fixed() +
  annotate('rect', xmin = ymin, xmax = ymax, ymin = xmin, ymax = xmax, fill = NA, colour = 'black', linewidth = 0.8) +
  guides(colour = guide_legend(override.aes = list(size = 3))) + theme_bw() +
  labs(x = 'y (horizontal in ImageDimPlot)', y = 'x (vertical)', colour = 'EBCC label',
       title = sprintf('mouse colon: dissection cut %g <= x <= %g, %g <= y <= %g', xmin, xmax, ymin, ymax),
       subtitle = sprintf('%d of %d cells; whole section is full-thickness colon wall, as is the human sample', sum(cells$cut), nrow(cells)))
dir.create('../figures/colon/spatial', showWarnings = FALSE, recursive = TRUE)
ggsave('../figures/colon/spatial/mouse_dissection.png', p, width = 12, height = 9, dpi = 120)
