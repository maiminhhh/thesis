# Human fetal liver, confident Seurat calls that EBCC relabels (same 7,211 cells as analysis/fl_seurat_overrides.R):
# one bar per change (Seurat label -> EBCC label), split by which of the other three annotators give EBCC's label.
# Changes with fewer than min_n cells pooled as 'Other changes'. Reads the obs export of analysis/fl_seurat_overrides.R.
# Usage: R_LIBS_SITE= R_LIBS_USER= <r45>/bin/Rscript --vanilla analysis/fl_seurat_overrides_flow.R
#   -> figures/cross_tissue/fl_seurat_overrides_flow.{png,pdf,csv}, _caption.txt
setwd(Sys.getenv('PROJECT_ROOT', '.'))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2) })

# ---- config ----
obs_csv = 'figures/cross_tissue/fl_seurat_overrides_obs.csv'
out = 'figures/cross_tissue/fl_seurat_overrides_flow'
others = c('SingleR' = 'celltype_singleR', 'Seurat (cross)' = 'celltype_seurat_cross', 'SingleR (cross)' = 'celltype_singleR_cross')
top_bins = 4:5
min_n = 30
expected_overrides = c(`4` = 3742, `5` = 3469)    # seurat_vs_ebcc/fetal_liver_human.csv, n_override

key = function(x) tolower(gsub('[^A-Za-z0-9]+', '_', x))
meta = read.csv(obs_csv, na.strings = 'NA')       # base parser (see analysis/fl_seurat_overrides.R)
set.seed(1)
held = !seq_len(nrow(meta)) %in% sample(nrow(meta), floor(nrow(meta) / 2))
m = meta[held, ]
keep = !is.na(m$celltype_seurat) & !is.na(m$labelscore_seurat)
m$bin = NA_integer_
m$bin[keep] = ntile(m$labelscore_seurat[keep], 5)
m$override = key(m$celltype_seurat) != key(m$ebcc_label)
chk = m |> filter(bin %in% top_bins) |> count(bin, wt = override)
stopifnot(all(chk$n == expected_overrides[as.character(chk$bin)]))

d = m |> filter(bin %in% top_bins, override)
agree = sapply(others, \(col) !is.na(d[[col]]) & key(d[[col]]) == key(d$ebcc_label))
d$backing = apply(agree, 1, \(x) if (all(x)) 'All three' else paste(names(others)[x], collapse = ' + '))
d$change = paste(d$celltype_seurat, '→', d$ebcc_label)
big = d |> count(change) |> filter(n >= min_n) |> pull(change)
d$change = if_else(d$change %in% big, d$change, 'Other changes')
main_backing = c('SingleR + SingleR (cross)', 'SingleR (cross)', 'Seurat (cross)', 'SingleR', 'All three')
d$backing = if_else(d$backing %in% main_backing, d$backing, 'Other pairs')
backing_lev = c(main_backing, 'Other pairs')
tab = d |> count(change, backing) |> group_by(change) |> mutate(total = sum(n)) |> ungroup()
write_csv(tab, paste0(out, '.csv'))
lev = tab |> distinct(change, total) |> arrange(change != 'Other changes', total) |> pull(change)
tot = tab |> distinct(change, total) |> mutate(pct = total / sum(total))

cols = c('SingleR + SingleR (cross)' = '#2A7F9E', 'SingleR (cross)' = '#E69F00', 'Seurat (cross)' = '#7B5EA7',
         'SingleR' = '#D2691E', 'All three' = '#1B3A57', 'Other pairs' = 'grey70')
p = ggplot(tab, aes(n, factor(change, levels = lev), fill = factor(backing, levels = backing_lev))) +
  geom_col(width = 0.72, colour = 'white', linewidth = 0.3) +
  geom_text(data = tot, aes(total, factor(change, levels = lev), label = sprintf('%s (%s)', format(total, big.mark = ','), if_else(pct < 0.01, '<1%', sprintf('%.0f%%', 100 * pct)))),
            inherit.aes = FALSE, hjust = -0.1, size = 3.3, colour = 'grey20', family = 'DejaVu Sans') +
  scale_fill_manual(values = cols, breaks = backing_lev, name = 'Annotators giving EBCC\'s label:') +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)), labels = scales::comma) +
  labs(x = 'Cells', y = 'Seurat label → EBCC label') +
  guides(fill = guide_legend(nrow = 2)) +
  theme_bw(base_size = 11, base_family = 'DejaVu Sans') +   # full-size arrow glyph
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = 'top',
        legend.justification = 'left')
ggsave(paste0(out, '.png'), p, width = 10, height = 5.5, dpi = 300)
ggsave(paste0(out, '.pdf'), p, width = 10, height = 5.5, device = cairo_pdf)
writeLines(paste(sprintf(paste('Human fetal liver, evaluation half: the %s cells in the top two fifths of Seurat\'s prediction score (>= %.3f)',
                               'whose EBCC label differs from Seurat\'s, by change of label; changes with fewer than %d cells pooled.'),
                         format(nrow(d), big.mark = ','), min(d$labelscore_seurat), min_n),
                 'Colour: which of the other three annotators (SingleR, Seurat cross-species, SingleR cross-species) give EBCC\'s label;',
                 'no cell lacks support from all three. A pruned SingleR call counts as not matching.', sep = '\n'), paste0(out, '_caption.txt'))
print(as.data.frame(tab |> arrange(desc(total), desc(n))))
