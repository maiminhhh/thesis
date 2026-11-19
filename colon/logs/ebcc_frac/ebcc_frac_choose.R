# Section 3b of ../EBCC_PIPELINE.md: choose frac (M fixed) by marker AUROC of the full EBCC fits.
# Marker score per cell and type = mean log-normalised expression of the type's panel markers
# (the 02 script's first `markers = list(...)`, as in ../analysis/visualisations.R); AUROC per type over the
# cells of one half; summary = mean over the types the labels call in that half.
# Cells split 50/50 (seed 1): frac with the best mean AUROC on the selection half (ties -> smaller
# frac), AUROC on the held-out half reported next to it.
# usage: Rscript ebcc_frac_choose.R human|mouse     (prints "CHANGED yes|no" as its last line)
which = commandArgs(TRUE)[1]
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
suppressMessages(source('../library.R'))   # Seurat, ggplot2 and the lab's library paths
file = if (which == 'human') '02_human_annotated.R' else '02_mouse_annotate.qmd'
src = readLines(file)
fracs = c(0, 0.25, 0.5, 1, 2, 5)
fits = lapply(setNames(fracs, fracs), function(f)
  readRDS(sprintf('data/02_annotated_data/ebcc_frac_fits/%s_frac%s.rds', which, format(f))))
M = unique(sapply(fits, `[[`, 'M')); stopifnot(length(M) == 1)

# same helpers as ../analysis/visualisations.R
read_markers = function(script) {
  l = readLines(script)
  s = grep('^markers = list', l)[1]
  e = s - 1 + which(l[s:length(l)] == ')')[1]
  env = new.env()
  eval(parse(text = l[s:e]), envir = env)
  env$markers
}
auroc = function(score, pos) {            # Mann-Whitney form, ties averaged
  r = rank(score); n1 = as.numeric(sum(pos)); n0 = as.numeric(sum(!pos))   # doubles: n1 * n0 overflows integers
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
key = function(x) tolower(gsub('[^A-Za-z0-9]+', '_', x))

obj = readRDS(sprintf('data/02_annotated_data/02_%s_annotated.rds', which))
rownames(obj) = obj[['Xenium']]@meta.data$gene_name
X = LayerData(obj, assay = 'Xenium', layer = 'data')
cells = colnames(obj); rm(obj); gc()
markers = lapply(read_markers(file), intersect, rownames(X))
markers = markers[lengths(markers) > 0]
score = sapply(markers, function(g) colMeans(as.matrix(X[g, , drop = FALSE])))   # cells x types
rm(X); gc()

set.seed(1)
half = rep('heldout', length(cells))
half[sample(length(cells), floor(length(cells) / 2))] = 'selection'

per_class = do.call(rbind, lapply(fracs, function(f) {
  fit = fits[[as.character(f)]]
  stopifnot(setequal(rownames(fit$posterior), cells))
  lab = key(as.character(fit$label[match(cells, rownames(fit$posterior))]))
  do.call(rbind, lapply(c('selection', 'heldout'), function(h) {
    in_h = half == h & !is.na(lab)
    called = unique(lab[in_h])
    do.call(rbind, lapply(names(markers), function(k) {
      if (!key(k) %in% called) return(NULL)
      pos = lab[in_h] == key(k)
      data.frame(frac = f, half = h, celltype = k, n_markers = length(markers[[k]]), n_cells = sum(pos),
                 auroc = if (all(pos)) NA_real_ else auroc(score[in_h, k], pos))
    }))
  }))
}))

summ = aggregate(auroc ~ frac + half, per_class, mean)
summ = reshape(summ, idvar = 'frac', timevar = 'half', direction = 'wide')
names(summ) = sub('^auroc\\.', 'mean_auroc_', names(summ))
summ$n_types = sapply(summ$frac, function(f) sum(per_class$frac == f & per_class$half == 'selection'))
summ = summ[order(summ$frac), c('frac', 'mean_auroc_selection', 'mean_auroc_heldout', 'n_types')]
best = max(summ$mean_auroc_selection)
pick = summ$frac[summ$mean_auroc_selection == best][1]   # fracs ascending: ties -> smaller frac
summ$chosen = summ$frac == pick
rownames(summ) = NULL
print(summ, digits = 4)

prev_M    = as.numeric(sub('^M = ([0-9.]+).*', '\\1', grep('^M = [0-9]', trimws(src), value = TRUE)))
prev_frac = as.numeric(sub('^frac = ([0-9.]+).*', '\\1', grep('^frac = [0-9]', trimws(src), value = TRUE)))
stopifnot(length(prev_M) == 1, length(prev_frac) == 1)
changed = M != prev_M || pick != prev_frac
rule = paste('M fixed from the held-out (M, frac) grid; full EBCC fit (R = 10) per frac;',
             'cells split 50/50 (set.seed(1)); frac with the highest mean marker AUROC on the selection half',
             '(ties -> smaller frac); frac tuned for marker agreement, so marker AUROC is not an independent check')
out = sprintf('data/02_annotated_data/%s_ebcc_frac_auroc', which)
write.csv(summ, paste0(out, '.csv'), row.names = FALSE)
write.csv(per_class, paste0(out, '_per_class.csv'), row.names = FALSE)
write.csv(data.frame(species = which, M = M, frac = pick, changed = changed,
                     previous_M = prev_M, previous_frac = prev_frac, rule = rule),
          paste0(out, '_chosen.csv'), row.names = FALSE)

d = rbind(data.frame(frac = summ$frac, half = 'selection', auroc = summ$mean_auroc_selection),
          data.frame(frac = summ$frac, half = 'held-out', auroc = summ$mean_auroc_heldout))
p = ggplot(d, aes(factor(frac), auroc, colour = half, group = half)) +
  geom_line() + geom_point(size = 2) +
  geom_vline(xintercept = match(pick, summ$frac), linetype = 'dashed') +
  labs(x = 'frac', y = 'mean marker AUROC (called types)', colour = NULL,
       title = sprintf('%s EBCC: marker AUROC vs frac (M = %g)', which, M),
       subtitle = sprintf('dashed: chosen frac = %g (best on the selection half)', pick)) +
  theme_minimal()
dir.create('../figures/colon/annotation', showWarnings = FALSE, recursive = TRUE)
ggsave(sprintf('../figures/colon/annotation/ebcc_frac_auroc_%s.png', which), p, width = 6, height = 4, dpi = 150)

if (pick == max(fracs)) message('NOTE: the largest frac won; the choice sits at the edge of the grid')
cat(sprintf('%s: chosen frac = %g with M = %g (selection AUROC %.4f, held-out %.4f; previously M = %g, frac = %g)\nCHANGED %s\n',
            which, pick, M, best, summ$mean_auroc_heldout[summ$chosen], prev_M, prev_frac, if (changed) 'yes' else 'no'))
