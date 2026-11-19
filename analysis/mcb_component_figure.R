# Appendix figure: Hsu's MCB (95% simultaneous) on the components of the overall score, 19 embedding methods (graph-based
# excluded, as section22<suffix>.pdf) x 4 tissues, outcome = improvement over unintegrated (dissected).
#   A batch correction   B bio-conservation   C bio-conservation without scANVI (18 methods; scANVI's row left blank)
# mcb() from analysis/mcb_analysis.R unchanged; rows, order, place numbers and input tiles as section22<suffix>.pdf (mean overall improvement).
# No footnotes in the image; caption in mcb_components_figure<suffix>_caption.txt.
# Usage: R_LIBS_SITE= R_LIBS_USER= <r45>/bin/Rscript --vanilla analysis/mcb_component_figure.R -> figures/mcb/mcb_components_figure<suffix>.{png,pdf}
setwd(Sys.getenv('PROJECT_ROOT', '.'))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr); library(ggplot2); library(patchwork) })

# ---- config ----
suffix = '_MAP_SAS'
graph_only = c('BBKNN', 'BBKNN UINMF', 'BBKNN scVI-latent', 'SAMap')   # as analysis/mcb_analysis.R
label_based = 'scANVI'
alpha = 0.05
input_levels = c('Native', 'scVI-latent', 'UINMF')
input_hl = c(Native = '#d9d9d9', `scVI-latent` = '#7B5EA7', UINMF = '#CCBB44')   # as analysis/mcb_analysis.R
input_text = c(Native = 'N', `scVI-latent` = 'S', UINMF = 'U')
out = sprintf('figures/mcb/mcb_components_figure%s', suffix)

ex = parse('analysis/mcb_analysis.R')
eval(ex[[which(sapply(ex, \(e) startsWith(deparse(e)[1], 'mcb = function')))]])
label_of = function(m) paste0(str_replace(m, '^(CCA|RPCA)', '\\1 (Seurat)'), ifelse(m == label_based, ' †', ''))
input_of = function(m) factor(case_when(str_detect(m, 'UINMF') ~ 'UINMF', str_detect(m, 'scVI-latent') ~ 'scVI-latent',
                                        TRUE ~ 'Native'), levels = input_levels)

d = read_csv(sprintf('figures/benchmark_summary/benchmark_improvement%s.csv', suffix), show_col_types = FALSE) |>
  filter(!method %in% .env$graph_only) |> select(method, tissue, overall = imp_overall, bio = imp_bio, batch = imp_batch)
stopifnot(n_distinct(d$method) == 19, nrow(d) == 76)
ord = mcb(d, 'overall')$res$method                                   # section 2.2 order: mean overall improvement
y_levels = rev(label_of(ord))

panels = list(A = list(data = d, y = 'batch', name = 'Batch correction'),
              B = list(data = d, y = 'bio', name = 'Bio-conservation'),
              C = list(data = filter(d, method != label_based), y = 'bio', name = 'Bio-conservation without scANVI'))
fits = lapply(panels, \(p) mcb(p$data, p$y))
summ = bind_rows(lapply(names(panels), function(k) {
  r = fits[[k]]; lead = r$res[which.max(r$res$mean), ]
  tibble(panel = k, score = panels[[k]]$name, n = nrow(r$res), could_be_best = sum(r$res$could_be_best), leader = lead$method,
         leader_mean = lead$mean, leader_lead = lead$D, d = r$d, se = r$se, threshold = r$d * r$se,
         leader_significant = lead$mcb_lo >= 0 & lead$D - r$d * r$se > 0)
}))

strip_df = tibble(method = ord, y = factor(label_of(ord), levels = y_levels), input = input_of(ord), pos = seq_along(ord))
strip = ggplot(strip_df, aes(1, y)) +
  geom_text(aes(x = 0, label = pos), size = 3.4, fontface = 'bold', colour = 'grey10') +
  geom_tile(aes(fill = input), width = 0.9, height = 0.85) +
  geom_text(aes(label = input_text[input], colour = input == 'Native'), size = 2.6, fontface = 'bold') +
  scale_x_continuous(limits = c(-0.5, 1.5)) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(values = input_hl, name = 'Input', drop = FALSE, labels = \(x) sprintf('%s (%s)', x, input_text[x])) +
  scale_colour_manual(values = c(`TRUE` = 'grey20', `FALSE` = 'white'), guide = 'none') +
  theme_void(base_size = 11) +
  theme(axis.text.y = element_text(colour = 'grey10', hjust = 1, size = 9.5, margin = margin(r = 2)))

mcb_panel = function(k) {
  r = fits[[k]]; s = summ[summ$panel == k, ]
  f = r$res |> mutate(y = factor(label_of(method), levels = y_levels))
  ggplot(f, aes(D, y)) +
    geom_vline(xintercept = 0, colour = 'grey30', linewidth = 0.5) +
    geom_errorbar(aes(xmin = mcb_lo, xmax = mcb_hi, colour = could_be_best), width = 0, orientation = 'y', linewidth = 0.9) +
    geom_point(aes(shape = could_be_best, colour = could_be_best), size = 2.6, fill = 'black', stroke = 0.9) +
    scale_y_discrete(drop = FALSE) +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), labels = c(`TRUE` = 'Could be best', `FALSE` = 'Worse than the best'),
                       name = NULL, breaks = c('TRUE', 'FALSE')) +
    scale_colour_manual(values = c(`TRUE` = 'black', `FALSE` = 'grey55'), labels = c(`TRUE` = 'Could be best',
                        `FALSE` = 'Worse than the best'), name = NULL, breaks = c('TRUE', 'FALSE')) +
    labs(x = 'Mean improvement minus\nthe best of the other methods', y = NULL,
         title = sprintf('%s\n(%d of %d could be best)', panels[[k]]$name, s$could_be_best, s$n)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = 'grey94'),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(face = 'bold', size = 11, hjust = 0.5, lineheight = 1.1))
}
p = strip + mcb_panel('A') + mcb_panel('B') + mcb_panel('C') +
  plot_layout(widths = c(0.16, 1, 1, 1), guides = 'collect') +
  plot_annotation(tag_levels = list(c('', 'A', 'B', 'C'))) &
  theme(legend.position = 'bottom', plot.tag = element_text(face = 'bold', size = 14))
ggsave(paste0(out, '.png'), p, width = 15, height = 8, dpi = 300, bg = 'white')
ggsave(paste0(out, '.pdf'), p, width = 15, height = 8, device = cairo_pdf, bg = 'white')
write_csv(summ, paste0(out, '_summary.csv'))
writeLines(paste(
  'Hsu\'s multiple comparisons with the best (95% simultaneous) for the two components of the overall score, 19 embedding-based methods',
  'x 4 tissues (graph-based methods excluded), outcome = improvement over unintegrated data. Point: a method\'s mean improvement minus',
  'the best mean of the other methods; interval constrained to include zero. Filled: could be the best (interval reaches above zero);',
  'open: significantly worse than the best. A: batch correction; B: bio-conservation; C: bio-conservation without scANVI (18 methods,',
  'its row left blank). Rows, order and numbers as in Figure section22 (mean overall improvement). † trained on the EBCC labels.',
  sprintf('Leader and margin (d x SE; a method could be best if it trails the leader by less): %s.',
          paste(sprintf('%s %s, %.3f (lead %.3f)', summ$panel, summ$leader, summ$threshold, summ$leader_lead), collapse = '; ')),
  sep = '\n'), paste0(out, '_caption.txt'))
print(as.data.frame(summ |> mutate(across(where(is.numeric), \(x) round(x, 4)))))
