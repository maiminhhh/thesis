# CCA, RPCA and FastMNN on native input as shared one-to-one orthologues are removed (dissected = 0%, then 20-80%).
# Scores loaded and computed exactly as analysis/orthologue_input.R (config + metric loading evaluated as written).
# Shown for the +MAP +SAS definition and the scib-only reference; improvement = overall minus Unintegrated, same tissue and case.
#   Rscript --vanilla analysis/native_shared_genes.R
# Output: figures/orthologue_input/native_shared_genes.{pdf,png,csv}
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr); library(stringr); library(ggplot2); library(patchwork)})
ex = parse('analysis/orthologue_input.R', keep.source = TRUE); hd = sapply(attr(ex, 'srcref'), \(s) as.character(s)[1])
stop_at = which(startsWith(hd, 'met = bind_rows(met, extra)')); stopifnot(length(stop_at) == 1)
for (i in seq_len(stop_at)) if (!grepl('^suppressPackageStartupMessages', hd[i])) eval(ex[[i]], globalenv())

show = c('CCA', 'RPCA', 'FastMNN'); method_cols = c(CCA = '#009988', RPCA = '#EE7733', FastMNN = '#7744AA')
sc = map_dfr(names(definitions), \(s) score_with(met, definitions[[s]]$bio, definitions[[s]]$batch) |>
               filter(method %in% c(show, reference)) |> mutate(definition = ifelse(s == '', 'scib only', '+MAP +SAS')))
stopifnot(nrow(sc) == 2 * 4 * length(cases) * 4)
ref = sc |> filter(method == reference) |> select(definition, tissue, case, ref = overall)
d = sc |> filter(method %in% show) |> left_join(ref, by = c('definition', 'tissue', 'case')) |>
  left_join(shared_genes, by = c('tissue', 'case')) |>
  mutate(removed = cases[case], improvement = overall - ref, method = factor(method, levels = show),
         tissue = factor(tissue, levels = names(tissue_labels))) |>
  select(definition, tissue, case, removed, genes, method, bio, batch, overall, unintegrated = ref, improvement)
write_csv(d, 'figures/orthologue_input/native_shared_genes.csv')

panel = function(y, ylab, def = '+MAP +SAS') {
  x = filter(d, definition == def)
  ggplot(x, aes(removed, .data[[y]], colour = method)) +
    { if (y == 'improvement') geom_hline(yintercept = 0, colour = 'grey55', linetype = 'dashed', linewidth = 0.4) } +
    geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
    facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(tissue_labels)) +
    scale_colour_manual(values = method_cols, name = 'Native input') +
    scale_x_continuous(breaks = cases, labels = paste0(cases, '%')) +
    labs(x = 'Shared one-to-one orthologues removed (0% = dissected)', y = ylab) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
                                     strip.text = element_text(face = 'bold', size = 10), legend.position = 'top')
}
p = (panel('overall', 'Overall score') + theme(axis.title.x = element_blank())) /
  panel('improvement', 'Improvement over\nUnintegrated (overall)') + plot_layout(guides = 'collect') &
  theme(legend.position = 'top')
ggsave('figures/orthologue_input/native_shared_genes.png', p, width = 10, height = 6, dpi = 300, bg = 'white')
ggsave('figures/orthologue_input/native_shared_genes.pdf', p, width = 10, height = 6, device = cairo_pdf, bg = 'white')

options(width = 200)
for (def in c('+MAP +SAS', 'scib only')) {
  cat('\n==== ', def, ': overall (improvement over Unintegrated) ====\n', sep = '')
  d |> filter(definition == def) |> mutate(v = sprintf('%.3f (%+.3f)', overall, improvement)) |>
    select(tissue, method, removed, v) |> pivot_wider(names_from = removed, values_from = v, names_prefix = 'rm') |>
    arrange(tissue, method) |> as.data.frame() |> print(row.names = FALSE)
  cat('Change in overall score, dissected -> 80% removed:\n')
  d |> filter(definition == def) |> group_by(tissue, method) |>
    summarise(change = overall[removed == 80] - overall[removed == 0], .groups = 'drop') |>
    pivot_wider(names_from = method, values_from = change) |> mutate(across(where(is.double), \(x) round(x, 3))) |>
    as.data.frame() |> print(row.names = FALSE)
}
cat('\nMean over tissues (+MAP +SAS):\n')
d |> filter(definition == '+MAP +SAS') |> group_by(method, removed) |> summarise(overall = mean(overall), impr = mean(improvement), .groups = 'drop') |>
  mutate(v = sprintf('%.3f (%+.3f)', overall, impr)) |> select(method, removed, v) |>
  pivot_wider(names_from = removed, values_from = v, names_prefix = 'rm') |> as.data.frame() |> print(row.names = FALSE)
