# Joint (M, frac) held-out grid for EBCC, used by 02_human_annotate.qmd and 02_mouse_annotate.qmd.
# heldout_ll(M, frac) is defined in each notebook: it fits EBCC on meta_ho (one annotator's label
# hidden for 10% of cells) and returns the log-probability of every hidden label.

GRID_RULE = paste('best = highest mean held-out log-likelihood; among the (M, frac) whose mean is',
                  'within 1 SE (sd/sqrt(n) of the best) of the best, take the smallest M, then the smallest frac')

# one forked process per frac, each fitting every M; `threads` BLAS/OpenMP threads per process
heldout_grid = function(heldout_ll, Ms = c(1, 2, 3, 5, 10), fracs = c(0, 0.25, 0.5, 1, 2), threads = 4) {
  per_frac = parallel::mclapply(fracs, function(f) {
    RhpcBLASctl::blas_set_num_threads(threads)
    RhpcBLASctl::omp_set_num_threads(threads)
    lapply(Ms, function(M) heldout_ll(M, f))
  }, mc.cores = length(fracs))
  bad = vapply(per_frac, function(x) inherits(x, 'try-error') || is.null(x), logical(1))
  if (any(bad)) stop('held-out fits failed for frac = ', paste(fracs[bad], collapse = ', '), ': ',
                     paste(sapply(per_frac[bad], as.character), collapse = '; '))

  grid = expand.grid(M = Ms, frac = fracs)          # M varies fastest, matching the unlist order below
  ll = unlist(per_frac, recursive = FALSE)
  grid$mean_ll   = vapply(ll, mean, numeric(1))
  grid$se        = vapply(ll, function(x) sd(x) / sqrt(length(x)), numeric(1))
  grid$n_heldout = lengths(ll)
  best = which.max(grid$mean_ll)
  grid$within_1se_of_best = grid$mean_ll >= grid$mean_ll[best] - grid$se[best]
  grid$best = seq_len(nrow(grid)) == best
  cand = grid[grid$within_1se_of_best, ]
  pick = cand[order(cand$M, cand$frac), ][1, ]
  grid$chosen = grid$M == pick$M & grid$frac == pick$frac
  list(ll = ll, grid = grid, M = pick$M, frac = pick$frac, rule = GRID_RULE)
}

# rds (all per-cell values + grid + rule), grid csv, chosen csv and heatmap; returns TRUE if the choice changed
save_heldout_grid = function(res, species, previous_M, previous_frac, ho = NULL,
                             data_dir = 'data/02_annotated_data', fig_dir = '../figures/breast_cancer/annotation') {
  base = file.path(data_dir, paste0(species, '_ebcc_heldout_grid'))
  saveRDS(c(res, list(ho = ho)), paste0(base, '.rds'))
  write.csv(res$grid, paste0(base, '.csv'), row.names = FALSE)
  changed = res$M != previous_M || res$frac != previous_frac
  write.csv(data.frame(species = species, M = res$M, frac = res$frac, changed = changed,
                       previous_M = previous_M, previous_frac = previous_frac),
            paste0(base, '_chosen.csv'), row.names = FALSE)

  g = res$grid
  g$M = factor(g$M); g$frac = factor(g$frac)
  p = ggplot(g, aes(frac, M, fill = mean_ll)) +
    geom_tile() +
    geom_text(aes(label = sprintf('%.4f', mean_ll)), size = 3) +
    geom_tile(data = g[g$within_1se_of_best, ], fill = NA, colour = 'black', linewidth = 0.4) +
    geom_tile(data = g[g$chosen, ], fill = NA, colour = 'black', linewidth = 1.6) +
    scale_fill_viridis_c(name = 'mean held-out\nlog-likelihood') +
    labs(title = sprintf('%s EBCC held-out grid: chosen M = %s, frac = %s', species, res$M, res$frac),
         subtitle = 'thin outline: within 1 SE of the best; thick outline: chosen') +
    theme_minimal()
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(fig_dir, paste0('ebcc_heldout_grid_', species, '.png')), p, width = 7, height = 5, bg = 'white')
  changed
}

# ---- M: the 1-SE rule applied to the frac = 0 column of the held-out grid ----
choose_M_frac0 = function(grid) {
  g = grid[grid$frac == 0, ]
  best = which.max(g$mean_ll)
  min(g$M[g$mean_ll >= g$mean_ll[best] - g$se[best]])
}

# ---- section 3b: frac by marker AUROC, M fixed (../EBCC_PIPELINE.md) ----
# markers: named list of marker genes (the 02 notebook's list); label_map: EBCC label -> marker group.
# Per cell and label: score = mean log-normalised expression of the label's markers; AUROC = P(a cell
# with that label scores higher than a cell without it); summary = mean over the labels called.
auroc = function(score, pos) {            # Mann-Whitney form, ties averaged
  r = rank(score); n1 = as.numeric(sum(pos)); n0 = as.numeric(sum(!pos))   # numeric: n1 * n0 overflows integers
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

frac_auroc = function(X, fits, markers, label_map, seed = 1) {
  markers = lapply(markers, intersect, rownames(X))
  label_map = label_map[label_map %in% names(markers)[lengths(markers) > 0]]
  scores = sapply(label_map, function(g) Matrix::colMeans(X[markers[[g]], , drop = FALSE]))
  set.seed(seed)
  sel = seq_len(ncol(X)) %in% sample(ncol(X), floor(ncol(X) / 2))
  halves = list(selection = sel, heldout = !sel)
  per_class = do.call(rbind, lapply(names(fits), function(f) {
    lab = as.character(fits[[f]]$label)
    do.call(rbind, lapply(names(halves), function(h) {
      idx = halves[[h]]
      do.call(rbind, lapply(names(label_map), function(k) {
        pos = lab[idx] == k
        data.frame(frac = as.numeric(f), half = h, celltype = k, marker_group = label_map[[k]],
                   n_cells = sum(pos),
                   auroc = if (sum(pos) == 0 || sum(!pos) == 0) NA_real_ else auroc(scores[idx, k], pos))
      }))
    }))
  }))
  summ = aggregate(auroc ~ frac + half, data = per_class, FUN = mean, na.action = na.omit)
  n_types = aggregate(auroc ~ frac + half, data = per_class, FUN = length, na.action = na.omit)
  wide = data.frame(frac = sort(unique(summ$frac)))
  wide$mean_auroc_selection = summ$auroc[summ$half == 'selection'][match(wide$frac, summ$frac[summ$half == 'selection'])]
  wide$mean_auroc_heldout   = summ$auroc[summ$half == 'heldout'][match(wide$frac, summ$frac[summ$half == 'heldout'])]
  wide$n_types = n_types$auroc[n_types$half == 'selection'][match(wide$frac, n_types$frac[n_types$half == 'selection'])]
  wide$chosen = seq_len(nrow(wide)) == which.max(wide$mean_auroc_selection)   # first max = smaller frac on ties
  list(per_class = per_class, summary = wide, frac = wide$frac[wide$chosen])
}

save_frac_auroc = function(res, species, M, previous_frac, previous_M,
                           data_dir = 'data/02_annotated_data', fig_dir = '../figures/breast_cancer/annotation') {
  base = file.path(data_dir, paste0(species, '_ebcc_frac_auroc'))
  write.csv(res$summary, paste0(base, '.csv'), row.names = FALSE)
  write.csv(res$per_class, paste0(base, '_per_class.csv'), row.names = FALSE)
  changed = M != previous_M || res$frac != previous_frac
  write.csv(data.frame(species = species, M = M, frac = res$frac, changed = changed,
                       previous_frac = previous_frac,
                       rule = paste('M: 1-SE rule on the frac = 0 column of the held-out grid;',
                                    'frac: highest mean marker AUROC on the selection half (50/50, seed 1), ties -> smaller frac;',
                                    'frac was tuned for marker agreement, so the marker AUROC is not an independent validation')),
            paste0(base, '_chosen.csv'), row.names = FALSE)
  long = rbind(data.frame(frac = res$summary$frac, half = 'selection', auroc = res$summary$mean_auroc_selection),
               data.frame(frac = res$summary$frac, half = 'held-out', auroc = res$summary$mean_auroc_heldout))
  p = ggplot(long, aes(factor(frac), auroc, colour = half, group = half)) +
    geom_line() + geom_point(size = 2) +
    geom_point(data = long[long$frac == res$frac, ], shape = 21, size = 5, colour = 'black', stroke = 1.2) +
    labs(x = 'frac', y = 'mean marker AUROC (called types)',
         title = sprintf('%s EBCC, M = %s: chosen frac = %s', species, M, res$frac),
         subtitle = 'circled: chosen on the selection half') +
    theme_minimal()
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(fig_dir, paste0('ebcc_frac_auroc_', species, '.png')), p, width = 7, height = 4.5, bg = 'white')
  changed
}
