# Dumbbell version of the ranking comparison with published scRNA-seq benchmarks (replaces literature_heatmap as Fig 4.1).
# Data and ranks are NOT recomputed here: the config + data/rank expressions of analysis/literature_comparison.R are evaluated
# as written (up to and including `pub`), so `ours` (rank / pct per tissue among the 11 basic methods) and `pub`
# (rank / pct per human-mouse task among shared methods; ties = half ranks; missing methods absent) are identical.
# Percentile p = (rank - 1) / (n - 1), 0 = best (analysis/literature_comparison.R's pct()).
#   Rscript --vanilla analysis/literature_dumbbell.R
# Output: figures/literature/literature_dumbbell.{pdf,png,csv}
suppressPackageStartupMessages({library(tidyverse); library(readxl)})

ex = parse('analysis/literature_comparison.R', keep.source = TRUE); hd = sapply(attr(ex, 'srcref'), \(s) as.character(s)[1])
stop_at = which(startsWith(hd, 'pub = pub_all'))
stopifnot(length(stop_at) == 1)
for (i in seq_len(stop_at)) if (!grepl('^source\\(', hd[i]) && !startsWith(hd[i], 'pub_all |>')) eval(ex[[i]], globalenv())

# ---------------- config ----------------
col_st  = '#C2410C'   # this study
col_pub = '#2563EB'   # published benchmarks
dy      = 0.12        # vertical offset of the faint per-tissue / per-task points
out_stem = file.path(out, 'literature_dumbbell')
task_names = c(`song:Pancreas_hs_mm` = 'BENGAL pancreas', `song:Heart_hs_mm` = 'BENGAL heart',
               `zhong:task4` = 'Zhong multi-tissue', `zhong:task21` = 'Zhong eye', `wang:task1` = 'Wang task 1')
# -----------------------------------------
lab = function(m) ifelse(m %in% c('CCA', 'RPCA'), paste0(m, ' (Seurat)'), m)   # as analysis/literature_comparison.R

stopifnot(length(basic) == 11, all(count(ours, tissue)$n == 11))
means = ours |> group_by(method) |> summarise(st = mean(pct)) |>
  left_join(pub |> group_by(method) |> summarise(pub = mean(pct), n_tasks = n()), by = 'method') |> arrange(st)
ord = means$method
yl  = setNames(rev(seq_along(ord)), ord)   # best (lowest ST mean) at the top

pts = bind_rows(ours |> transmute(method, pct, grp = 'st'), pub |> transmute(method, pct, grp = 'pub')) |>
  mutate(y = yl[method] + ifelse(grp == 'st', dy, -dy), key = paste(grp, 'each'))
mpts = means |> pivot_longer(c(st, pub), names_to = 'grp', values_to = 'pct') |> filter(!is.na(pct)) |>
  mutate(y = yl[method], key = paste(grp, 'mean'))

keys = c(`st each` = 'This study (ST): per tissue', `st mean` = 'This study (ST): mean',
         `pub each` = 'scRNA-seq benchmarks: per task', `pub mean` = 'scRNA-seq benchmarks: mean')
kcol  = c(`st each` = alpha(col_st, 0.35), `st mean` = col_st, `pub each` = alpha(col_pub, 0.35), `pub mean` = col_pub)
ksize = c(`st each` = 2.1, `st mean` = 3.4, `pub each` = 2.1, `pub mean` = 3.4)

fig = ggplot() +
  geom_segment(data = filter(means, !is.na(pub)), aes(x = st, xend = pub, y = yl[method], yend = yl[method]),
               colour = 'grey60', linewidth = 0.7) +
  geom_point(data = pts, aes(pct, y, colour = key, size = key), stroke = 0) +
  geom_point(data = mpts, aes(pct, y, colour = key, size = key)) +
  scale_colour_manual(values = kcol, breaks = names(keys), labels = keys, name = NULL) +
  scale_size_manual(values = ksize, breaks = names(keys), labels = keys, name = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), size = guide_legend(nrow = 2, byrow = TRUE)) +
  scale_y_continuous(breaks = yl, labels = lab(names(yl)), expand = expansion(add = 0.5)) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c('0 (best)', '0.5', '1 (worst)'),
                     expand = expansion(add = 0.03)) +
  labs(x = 'Rank percentile', y = NULL) +
  theme_minimal(base_size = 9.5) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey90', linewidth = 0.4),
        axis.text = element_text(colour = 'grey10', size = 9), axis.title.x = element_text(size = 9.5),
        legend.position = 'top', legend.text = element_text(size = 8.5), legend.key.width = unit(10, 'pt'),
        legend.margin = margin(0, 0, 0, 0), plot.margin = margin(4, 8, 4, 4))
ggsave(paste0(out_stem, '.pdf'), fig, width = 6.3, height = 4.5, device = cairo_pdf, bg = 'white')
ggsave(paste0(out_stem, '.png'), fig, width = 6.3, height = 4.5, dpi = 300, bg = 'white')

# ---- check table: every percentile and both means ----
wide = ours |> transmute(method, col = tissue_labels[tissue], pct) |>
  bind_rows(pub |> transmute(method, col = task_names[paste(paper, task, sep = ':')], pct)) |>
  pivot_wider(names_from = col, values_from = pct) |>
  select(method, all_of(unname(tissue_labels)), all_of(unname(task_names)))
tab = means |> select(method, ST_mean = st, pub_mean = pub, n_tasks) |> left_join(wide, by = 'method')
write_csv(tab, paste0(out_stem, '.csv'))
cat('Rank percentiles (0 = best), rows in plot order; NA = method not in that benchmark task\n')
print(as.data.frame(tab |> mutate(across(where(is.double), \(x) round(x, 3)))), row.names = FALSE)
cat('n per task:\n'); print(as.data.frame(distinct(pub, paper, task, n)), row.names = FALSE)
cat('wrote', paste0(out_stem, '.{pdf,png,csv}\n'))
