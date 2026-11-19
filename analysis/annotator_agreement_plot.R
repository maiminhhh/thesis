# Agreement heatmap between two annotators of one dataset, drawn with plot_annotation_agreement() from analysis/visualisations.R
# (as seurat_singleR_agreement.png). Counts from analysis/annotator_agreement_export.py.
#   R_LIBS_SITE= R_LIBS_USER= .../envs/r45/bin/Rscript --vanilla analysis/annotator_agreement_plot.R <in.csv> <out stem> "<title>" "<x title>" "<y title>"
suppressPackageStartupMessages({library(tidyverse); library(patchwork)})
a = commandArgs(TRUE)
vis = parse('analysis/visualisations.R'); vis_h = sapply(vis, \(e) deparse(e)[1])
for (f in c('save_plot = function', 'plot_annotation_agreement = function')) eval(vis[[which(startsWith(vis_h, f))]], globalenv())
res = read_csv(a[1], show_col_types = FALSE)
p = plot_annotation_agreement(res, title = a[3], axis_titles = c(a[4], a[5]))
save_plot(p, paste0(a[2], '.png'), width = 11, height = 8)
ggsave(paste0(a[2], '.pdf'), p, width = 11, height = 8, device = cairo_pdf)
cat('wrote', paste0(a[2], '.{png,pdf}\n'))
