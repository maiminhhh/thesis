# Do Seurat CCA / RPCA lose improvement more slowly than other algorithms on native input as shared orthologues are removed?
# Scores loaded exactly as analysis/orthologue_input.R (config + metric loading evaluated as written); BBKNN scored on its own metrics.
# impr = overall minus Unintegrated (same tissue, case). Model per definition, mirroring the 3.2.4 genes model
# (input -> algorithm, the factor of interest, so fixed; algorithm*tissue as input*tissue there):
#   impr ~ algorithm * n + algorithm * tissue + (1 | dataset); n = shared one-to-one orthologues in the dataset (raw count),
#   tissue fixed (sum-to-zero), so the slope on n is estimated within tissues; dataset = tissue x case.
# Per-algorithm slopes per shared gene (emtrends, Kenward-Roger; larger = loses more improvement as genes are removed),
# Seurat (CCA, RPCA) vs the other algorithms, and Seurat vs each other algorithm (unadjusted).
# Algorithm sets: section 2.4 (Harmony, BBKNN, FastMNN, RPCA, CCA, Scanorama) and all 10 native algorithms.
# Also: share of the dissected improvement kept at 80% removed (descriptive, per tissue).
#   Rscript --vanilla analysis/native_slopes.R
# Output: figures/orthologue_input/native_slopes<tag>.{csv,pdf,png}, native_slopes_contrasts<tag>.csv (tag '_no80' with drop80); printout.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr); library(stringr); library(ggplot2); library(patchwork); library(emmeans)})
ex = parse('analysis/orthologue_input.R', keep.source = TRUE); hd = sapply(attr(ex, 'srcref'), \(s) as.character(s)[1])
stop_at = which(startsWith(hd, 'met = bind_rows(met, extra)')); stopifnot(length(stop_at) == 1)
for (i in seq_len(stop_at)) if (!grepl('^suppressPackageStartupMessages', hd[i])) eval(ex[[i]], globalenv())
emm_options(lmer.df = 'kenward-roger', lmerTest.limit = 1e5, pbkrtest.limit = 1e5)
options(width = 220, contrasts = c('contr.sum', 'contr.poly'))

sets = list(`section 2.4` = c('CCA', 'RPCA', 'Harmony', 'FastMNN', 'Scanorama', 'BBKNN'),
            `all native` = c('CCA', 'RPCA', 'Harmony', 'FastMNN', 'Scanorama', 'BBKNN', 'LIGER', 'SAMap', 'scVI', 'scANVI'))
seurat = c('CCA', 'RPCA')
drop80 = 'drop80' %in% commandArgs(TRUE)   # Rscript analysis/native_slopes.R drop80: refit without the 80%-removed case (genes8)
tag = if (drop80) '_no80' else ''
slopes = list(); cons = list(); kept = list()
for (s in names(definitions)) {
  def = definitions[[s]]; dl = ifelse(s == '', 'scib only', '+MAP +SAS')
  sc = score_with(met, def$bio, def$batch)
  ref = sc |> filter(method == reference) |> select(tissue, case, ref = overall)
  d0 = sc |> filter(method %in% sets$`all native`) |> left_join(ref, by = c('tissue', 'case')) |>
    left_join(shared_genes, by = c('tissue', 'case')) |>
    mutate(impr = overall - ref, n = genes, r10 = cases[case] / 10, dataset = paste(tissue, case), tissue = factor(tissue))
  stopifnot(nrow(d0) == 10 * 20, !anyNA(d0$impr))
  kept_d0 = d0
  if (drop80) d0 = filter(d0, case != 'genes8')
  kept[[dl]] = kept_d0 |> group_by(method, tissue) |> summarise(i0 = impr[r10 == 0], i80 = impr[r10 == 8], .groups = 'drop') |>
    mutate(definition = dl)
  for (st in names(sets)) {
    d = d0 |> filter(method %in% sets[[st]]) |> mutate(algorithm = factor(method, levels = sets[[st]]))
    f = suppressMessages(lmerTest::lmer(impr ~ algorithm * n + algorithm * tissue + (1 | dataset), data = d))
    tr = emtrends(f, ~ algorithm, var = 'n')
    slopes[[length(slopes) + 1]] = as.data.frame(summary(tr, infer = TRUE)) |> mutate(definition = dl, set = st,
      at123 = as.data.frame(emmeans(f, ~ algorithm, at = list(n = 123)))$emmean,
      at25 = as.data.frame(emmeans(f, ~ algorithm, at = list(n = 25)))$emmean)
    lv = levels(d$algorithm); w = function(a) as.numeric(lv %in% a) / length(a)
    L = c(list(`Seurat (CCA, RPCA) - others` = w(seurat) - w(setdiff(lv, seurat))),
          setNames(lapply(setdiff(lv, seurat), \(o) w(seurat) - w(o)), paste('Seurat -', setdiff(lv, seurat))),
          setNames(lapply(setdiff(lv, seurat), \(o) w('CCA') - w(o)), paste('CCA -', setdiff(lv, seurat))),
          setNames(lapply(setdiff(lv, seurat), \(o) w('RPCA') - w(o)), paste('RPCA -', setdiff(lv, seurat))))
    cons[[length(cons) + 1]] = as.data.frame(summary(contrast(tr, L), infer = TRUE)) |> mutate(definition = dl, set = st)
  }
}
slopes = bind_rows(slopes) |> select(definition, set, algorithm, at123, at25, slope = n.trend, SE, df, lower.CL, upper.CL)
cons = bind_rows(cons) |> select(definition, set, contrast, estimate, SE, df, lower.CL, upper.CL, p.value)
write_csv(slopes, sprintf('figures/orthologue_input/native_slopes%s.csv', tag)); write_csv(cons, sprintf('figures/orthologue_input/native_slopes_contrasts%s.csv', tag))

r = \(x) round(x, 4)
for (dl in c('+MAP +SAS', 'scib only')) for (st in names(sets)) {
  cat(sprintf('\n======== %s | %s ========\nSlope = change in improvement over Unintegrated per shared gene (larger = loses more as genes are removed); fitted improvement at 123 and 25 shared genes (mean over tissues)\n', dl, st))
  slopes |> filter(definition == dl, set == st) |> arrange(desc(slope)) |> mutate(across(where(is.double), r)) |>
    select(-definition, -set, -df) |> as.data.frame() |> print(row.names = FALSE)
  cat('Contrasts of slopes per shared gene (negative = Seurat loses improvement more slowly):\n')
  cons |> filter(definition == dl, set == st, str_starts(contrast, 'Seurat')) |> mutate(across(where(is.double), r)) |>
    select(-definition, -set) |> as.data.frame() |> print(row.names = FALSE)
}
cat('\n======== Observed improvement over Unintegrated, dissected -> 80% removed, per tissue ========\n')
bind_rows(kept) |> mutate(v = sprintf('%+.3f -> %+.3f', i0, i80)) |> select(definition, method, tissue, v) |>
  pivot_wider(names_from = tissue, values_from = v) |> arrange(definition, method) |> as.data.frame() |> print(row.names = FALSE)

# figure: slopes with 95% CI, both definitions, all native algorithms
pd = slopes |> filter(set == 'all native') |> group_by(definition) |> mutate(algorithm = as.character(algorithm)) |> ungroup()
ord = pd |> filter(definition == '+MAP +SAS') |> arrange(desc(slope)) |> pull(algorithm)
pd = pd |> mutate(algorithm = factor(algorithm, levels = ord), seurat = algorithm %in% seurat)
p = ggplot(pd, aes(slope, algorithm, colour = seurat)) +
  geom_vline(xintercept = 0, colour = 'grey55', linetype = 'dashed', linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.25, linewidth = 0.6, orientation = 'y') + geom_point(size = 2.2) +
  facet_wrap(~ definition) +
  scale_colour_manual(values = c(`TRUE` = '#C2410C', `FALSE` = 'grey30'), labels = c(`TRUE` = 'Seurat (CCA, RPCA)', `FALSE` = 'Other'), name = NULL) +
  labs(x = 'Slope of improvement over Unintegrated per shared orthologue (95% CI)\nlarger = loses more improvement as shared genes are removed', y = NULL) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
                                   strip.background = element_blank(), strip.text = element_text(face = 'bold'), legend.position = 'top')
ggsave(sprintf('figures/orthologue_input/native_slopes%s.png', tag), p, width = 8, height = 4.2, dpi = 300, bg = 'white')
ggsave(sprintf('figures/orthologue_input/native_slopes%s.pdf', tag), p, width = 8, height = 4.2, device = cairo_pdf, bg = 'white')
