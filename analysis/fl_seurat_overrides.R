# Human fetal liver: where Seurat (the best annotator) is confident but EBCC relabels the cell, what backs EBCC's label?
# Cells: evaluation half (set.seed(1) split, as in analysis/vote_vs_ebcc.R / analysis/seurat_vs_ebcc.R), top two fifths of Seurat's
# prediction score (dplyr::ntile over the cells Seurat labels, as in analysis/seurat_vs_ebcc.R), EBCC label != Seurat label.
#   A. how many of the other three annotators (SingleR, Seurat cross, SingleR cross) give EBCC's label (0-3;
#      a pruned SingleR call counts as not matching)
#   B. EBCC's new labels, stacked by A
# Labels read from obs of the saved 02 h5ad (same cells, same order as the .rds; checked against seurat_vs_ebcc's
# override counts below), because Seurat no longer loads under R 4.5. Nothing is refitted or rescored.
# Usage: python3 export below, then R_LIBS_SITE= R_LIBS_USER= <r45>/bin/Rscript --vanilla analysis/fl_seurat_overrides.R
#   -> figures/cross_tissue/fl_seurat_overrides.{png,pdf}, _counts.csv, _labels.csv, _caption.txt
setwd(Sys.getenv('PROJECT_ROOT', '.'))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork) })

# ---- config ----
h5ad = 'fetal_liver/data/02_annotated_data/02_human_annotated.h5ad'
obs_csv = 'figures/cross_tissue/fl_seurat_overrides_obs.csv'
out = 'figures/cross_tissue/fl_seurat_overrides'
others = c('SingleR' = 'celltype_singleR', 'Seurat (cross)' = 'celltype_seurat_cross', 'SingleR (cross)' = 'celltype_singleR_cross')
top_bins = 4:5
expected_overrides = c(`4` = 3742, `5` = 3469)    # seurat_vs_ebcc/fetal_liver_human.csv, n_override

if (!file.exists(obs_csv)) system2('python3', c('-c', shQuote(sprintf(paste(
  "import h5py,numpy as np,csv; o=h5py.File('%s','r')['obs']",
  "cols=['celltype_seurat','celltype_singleR','celltype_seurat_cross','celltype_singleR_cross','ebcc_label','labelscore_seurat']",
  "g=lambda d: np.where(d['codes'][:]<0,'NA',d['categories'][:].astype(str)[np.maximum(d['codes'][:],0)]) if isinstance(d,h5py.Group) else d[:]",
  "D=[g(o[c]) for c in cols]; w=csv.writer(open('%s','w',newline='')); w.writerow(cols); w.writerows(zip(*D))", sep = '\n'), h5ad, obs_csv))))
key = function(x) tolower(gsub('[^A-Za-z0-9]+', '_', x))
meta = read.csv(obs_csv, na.strings = 'NA')   # base parser: readr's is off by an ulp near 1, which moves cells across the 4/5 cut

set.seed(1)
held = !seq_len(nrow(meta)) %in% sample(nrow(meta), floor(nrow(meta) / 2))
m = meta[held, ]
keep = !is.na(m$celltype_seurat) & !is.na(m$labelscore_seurat)
m$bin = NA_integer_
m$bin[keep] = ntile(m$labelscore_seurat[keep], 5)
m$override = key(m$celltype_seurat) != key(m$ebcc_label)
chk = m |> filter(bin %in% top_bins) |> count(bin, wt = override)
stopifnot(all(chk$n == expected_overrides[as.character(chk$bin)]))   # same cells as analysis/seurat_vs_ebcc.R

d = m |> filter(bin %in% top_bins, override)
d$support = rowSums(sapply(others, \(col) !is.na(d[[col]]) & key(d[[col]]) == key(d$ebcc_label)))
for (nm in names(others)) d[[nm]] = !is.na(d[[others[[nm]]]]) & key(d[[others[[nm]]]]) == key(d$ebcc_label)

counts = d |> count(support) |> complete(support = 0:3, fill = list(n = 0)) |> mutate(pct = n / sum(n))
labels = d |> count(ebcc_label, support) |> group_by(ebcc_label) |> mutate(total = sum(n)) |> ungroup() |>
  arrange(desc(total), support)
transitions = d |> count(celltype_seurat, ebcc_label, sort = TRUE) |> mutate(pct = n / sum(n))
per_annotator = sapply(names(others), \(nm) sum(d[[nm]]))
write_csv(counts, paste0(out, '_counts.csv'))
write_csv(labels, paste0(out, '_labels.csv'))
write_csv(transitions, paste0(out, '_transitions.csv'))

sup_cols = c('0' = '#D9D9D9', '1' = '#9ECAE1', '2' = '#4292C6', '3' = '#08519C')
pA = ggplot(counts, aes(factor(support), n, fill = factor(support))) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf('%s\n(%.0f%%)', format(n, big.mark = ','), 100 * pct)), vjust = -0.25, size = 3.2,
            colour = 'grey20', lineheight = 0.9) +
  scale_fill_manual(values = sup_cols, guide = 'none') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::comma) +
  labs(x = 'Other annotators giving EBCC\'s label (of SingleR, Seurat cross, SingleR cross)', y = 'Cells')
lev = labels |> distinct(ebcc_label, total) |> arrange(total) |> pull(ebcc_label)
tot = labels |> distinct(ebcc_label, total)
pB = ggplot(labels, aes(n, factor(ebcc_label, levels = lev), fill = factor(support, levels = 3:0))) +
  geom_col(width = 0.7, colour = 'white', linewidth = 0.25) +
  geom_text(data = tot, aes(total, factor(ebcc_label, levels = lev), label = format(total, big.mark = ',')),
            inherit.aes = FALSE, hjust = -0.15, size = 3.2, colour = 'grey20') +
  scale_fill_manual(values = sup_cols, breaks = as.character(0:3), name = 'Other annotators agreeing') +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12)), labels = scales::comma) +
  labs(x = 'Cells', y = 'EBCC\'s label')
p = (pA | pB) + plot_layout(widths = c(1, 1.3)) + plot_annotation(tag_levels = 'A') &
  theme_bw(base_size = 11) &
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), legend.position = 'top',
        plot.tag = element_text(face = 'bold', size = 14))
p[[2]] = p[[2]] + theme(panel.grid.major.x = element_line(colour = 'grey92'), panel.grid.major.y = element_blank())
ggsave(paste0(out, '.png'), p, width = 12, height = 5, dpi = 300)
ggsave(paste0(out, '.pdf'), p, width = 12, height = 5, device = cairo_pdf)
writeLines(paste(sprintf(paste('Human fetal liver, evaluation half: the %s cells in the top two fifths of Seurat\'s prediction score',
                               '(score >= %.3f) whose EBCC label differs from Seurat\'s.'), format(nrow(d), big.mark = ','),
                         min(d$labelscore_seurat)),
                 'A: number of the other three annotators (SingleR, Seurat cross-species, SingleR cross-species) whose label is EBCC\'s label;',
                 'a pruned SingleR call counts as not matching. B: EBCC\'s labels for these cells, coloured by the same count; totals at the bar ends.',
                 sep = '\n'), paste0(out, '_caption.txt'))

print(as.data.frame(counts)); print(per_annotator)
print(as.data.frame(labels |> select(-total) |> pivot_wider(names_from = support, values_from = n, values_fill = 0) |>
  mutate(total = rowSums(across(-ebcc_label))) |> arrange(desc(total))))
print(head(as.data.frame(transitions), 12))
cat('Seurat score range:', range(d$labelscore_seurat), '\n')
