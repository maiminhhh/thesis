# EBCC model selection on a joint (M, frac) grid, shared by 02_human_annotate.R and
# 02_mouse_annotate.R. Each (M, frac) fit's held-out log-likelihoods are cached in
# data/02_annotated_data/ebcc_grid/<species>/, so the grid can be split across processes:
#   EBCC_FRACS=0.25,0.5  compute only these frac values, then stop (worker mode)
# The cache is cleared when the annotations, held-out set or trusted classes change.

ebcc_grid_dir = function(species) file.path('data/02_annotated_data/ebcc_grid', species)

# log P(hidden label) for every held-out label, from one fit on the masked data
ebcc_heldout_one = function(M, frac, meta_ho, ho, ann_cols, make_beta, R = 3) {
  f = run_ebcc(meta_ho, ann_cols, M = M, R = R, keep_rho = TRUE,
               beta = make_beta(M, frac), n_cores = R)
  unlist(lapply(ann_cols, function(a) {
    s = ho[ho$ann == a, ]
    l = match(s$truth, f$classes); ok = !is.na(l)
    # P(hidden label) = sum over class-subtypes of P(cell in it) x P(annotator gives that label)
    log(rowSums(f$rho[s$cell[ok], , drop = FALSE] *
                t(f$par$v[[a]][, l[ok], drop = FALSE])))
  }))
}

# fit every (M, frac) not yet cached; in worker mode only the EBCC_FRACS values
run_ebcc_grid = function(species, Ms, fracs, meta_ho, ho, ann_cols, make_beta, trusted) {
  dir = ebcc_grid_dir(species)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  key = rlang::hash(list(meta_ho[ann_cols], ho, trusted))
  key_file = file.path(dir, 'key.txt')
  if (!file.exists(key_file) || readLines(key_file) != key) {
    unlink(list.files(dir, pattern = '^M.*\\.rds$', full.names = TRUE))
    writeLines(key, key_file)
  }
  todo = expand_grid(M = Ms, frac = fracs)
  only = Sys.getenv('EBCC_FRACS')
  if (nzchar(only)) todo = filter(todo, frac %in% as.numeric(strsplit(only, ',')[[1]]))
  for (i in seq_len(nrow(todo))) {
    f = file.path(dir, sprintf('M%s_frac%s.rds', todo$M[i], todo$frac[i]))
    if (file.exists(f)) next
    message(sprintf('%s: held-out fit M = %s, frac = %s', species, todo$M[i], todo$frac[i]))
    saveRDS(ebcc_heldout_one(todo$M[i], todo$frac[i], meta_ho, ho, ann_cols, make_beta), f)
  }
  if (nzchar(only)) quit(save = 'no')   # worker mode: selection runs once all fracs are cached
}

# rule: among (M, frac) whose mean held-out ll is within 1 SE (of the best's mean) of
# the best, the smallest M, then the smallest frac
select_ebcc_grid = function(species, Ms, fracs, previous, fig_dir = '../figures/fetal_liver/annotation') {
  dir = ebcc_grid_dir(species)
  grid = expand_grid(M = Ms, frac = fracs)
  values = purrr::map2(grid$M, grid$frac, function(M, frac) {
    f = file.path(dir, sprintf('M%s_frac%s.rds', M, frac))
    if (!file.exists(f)) stop('missing held-out fit ', f)
    readRDS(f)
  })
  grid$mean_ll   = purrr::map_dbl(values, mean)
  grid$se        = purrr::map_dbl(values, \(v) sd(v) / sqrt(length(v)))
  grid$n_heldout = lengths(values)
  best = which.max(grid$mean_ll)
  grid$within_1se_of_best = grid$mean_ll >= grid$mean_ll[best] - grid$se[best]
  grid$best = seq_len(nrow(grid)) == best
  chosen = grid |> mutate(i = row_number()) |> filter(within_1se_of_best) |>
    arrange(M, frac) |> slice(1) |> pull(i)
  grid$chosen = seq_len(nrow(grid)) == chosen
  M = grid$M[chosen]; frac = grid$frac[chosen]
  changed = M != previous$M || frac != previous$frac
  rule = 'within 1 SE of the best mean held-out ll (SE of the best); smallest M, then smallest frac'
  print(as.data.frame(grid))
  message(sprintf('%s: best M = %s, frac = %s; chosen M = %s, frac = %s (%s)',
                  species, grid$M[best], grid$frac[best], M, frac,
                  if (changed) sprintf('changed from M = %s, frac = %s', previous$M, previous$frac)
                  else 'unchanged'))

  out = file.path('data/02_annotated_data', paste0(species, '_ebcc_heldout_grid'))
  saveRDS(list(values = setNames(values, sprintf('M%s_frac%s', grid$M, grid$frac)),
               summary = grid, rule = rule), paste0(out, '.rds'))
  write_csv(grid, paste0(out, '.csv'))
  write_csv(tibble(species = species, M = M, frac = frac, changed = changed,
                   previous_M = previous$M, previous_frac = previous$frac),
            paste0(out, '_chosen.csv'))

  p = ggplot(grid, aes(factor(frac), factor(M), fill = mean_ll)) +
    geom_tile(colour = 'white') +
    geom_tile(data = filter(grid, within_1se_of_best), fill = NA, colour = 'black', linewidth = 0.4) +
    geom_tile(data = filter(grid, chosen), fill = NA, colour = 'black', linewidth = 1.6) +
    geom_text(aes(label = sprintf('%.4f', mean_ll)), size = 3.5) +
    scale_fill_viridis_c(name = 'Mean held-out\nlog-likelihood') +
    labs(title = sprintf('%s: EBCC held-out grid', species),
         subtitle = sprintf('thin outline: within 1 SE of the best; thick: chosen (M = %s, frac = %s)', M, frac),
         x = 'frac (trusted-class prior strength)', y = 'M (subtypes per class)') +
    theme_classic(base_size = 13)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(fig_dir, sprintf('ebcc_heldout_grid_%s.png', species)), p,
         width = 8, height = 6, dpi = 300, bg = 'white')
  list(M = M, frac = frac, changed = changed)
}

# ---- frac by marker AUROC, M fixed (EBCC_PIPELINE.md section 3b) ----
# Full fits (R = 10, all cells) for each frac are cached in
# data/02_annotated_data/ebcc_frac/<species>/, so they can run as separate processes:
#   EBCC_FRAC_FITS=0.5,1  fit only these frac values, then stop (worker mode)
#   EBCC_NCORES=5         parallel restarts per fit (default 10)
ebcc_frac_dir = function(species) file.path('data/02_annotated_data/ebcc_frac', species)

run_ebcc_frac_fits = function(species, meta, ann_cols, M, fracs, make_beta, trusted) {
  dir = ebcc_frac_dir(species)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  key = rlang::hash(list(meta[ann_cols], M, trusted))
  key_file = file.path(dir, 'key.txt')
  if (!file.exists(key_file) || readLines(key_file) != key) {
    unlink(list.files(dir, pattern = '^frac.*\\.rds$', full.names = TRUE))
    writeLines(key, key_file)
  }
  todo = fracs
  only = Sys.getenv('EBCC_FRAC_FITS')
  if (nzchar(only)) todo = intersect(fracs, as.numeric(strsplit(only, ',')[[1]]))
  n_cores = as.integer(Sys.getenv('EBCC_NCORES', '10'))
  for (frac in todo) {
    f = file.path(dir, sprintf('frac%s.rds', frac))
    if (file.exists(f)) next
    message(sprintf('%s: full EBCC fit M = %s, frac = %s', species, M, frac))
    fit = run_ebcc(meta, ann_cols, M = M, R = 10, beta = make_beta(M, frac), n_cores = n_cores)
    saveRDS(list(label = fit$label, posterior = fit$posterior, par = fit$par,
                 elbo_restarts = fit$elbo_restarts, M = M, frac = frac), f)
  }
  if (nzchar(only)) quit(save = 'no')   # worker mode: selection runs once all fits exist
}

# P(a cell labelled with the type has a higher marker score than one labelled otherwise),
# Mann-Whitney form with ties averaged
auroc = function(score, pos) {
  n1 = as.numeric(sum(pos)); n0 = as.numeric(sum(!pos))   # numeric: n1 * n0 overflows integer
  (sum(rank(score)[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# per cell type and half: AUROC of the type's mean marker expression (log-normalised)
# for cells labelled with the type vs all others; only types the labels call
marker_auroc_halves = function(labels, scores, half) {
  purrr::map_dfr(c('selection', 'heldout'), function(h) {
    keep = half == h
    lab = labels[keep]
    purrr::map_dfr(colnames(scores), function(k) {
      pos = lab == k
      if (!any(pos) || all(pos)) return(NULL)
      tibble(half = h, celltype = k, n_cells = sum(pos), auroc = auroc(scores[keep, k], pos))
    })
  })
}

select_frac_by_auroc = function(species, obj, M, fracs, markers, previous,
                                fig_dir = '../figures/fetal_liver/annotation') {
  fits = purrr::map(fracs, function(frac) {
    f = file.path(ebcc_frac_dir(species), sprintf('frac%s.rds', frac))
    if (!file.exists(f)) stop('missing EBCC fit ', f)
    readRDS(f)
  })
  # marker score per cell and type: mean log-normalised expression of the type's panel markers
  X = LayerData(obj, assay = 'Xenium', layer = 'data')
  markers = lapply(markers, intersect, rownames(X))
  markers = markers[lengths(markers) > 0]
  scores = sapply(markers, function(g) Matrix::colMeans(X[g, , drop = FALSE]))
  rownames(scores) = colnames(X)
  set.seed(1)
  half = rep('heldout', ncol(X))
  half[sample(ncol(X), floor(ncol(X) / 2))] = 'selection'

  per_class = purrr::map2_dfr(fits, fracs, function(fit, frac) {
    lab = as.character(fit$label)[match(colnames(X), names(fit$label) %||% rownames(fit$posterior))]
    marker_auroc_halves(lab, scores, half) |> mutate(frac = frac, .before = 1)
  })
  summary = per_class |>
    group_by(frac, half) |>
    summarise(mean_auroc = mean(auroc), n_types = n(), .groups = 'drop') |>
    pivot_wider(names_from = half, values_from = c(mean_auroc, n_types)) |>
    arrange(frac)
  # highest mean AUROC on the selection half; ties -> smaller frac (fracs are ascending)
  chosen_i = which.max(round(summary$mean_auroc_selection, 12))
  summary$chosen = seq_len(nrow(summary)) == chosen_i
  frac = summary$frac[chosen_i]
  changed = M != previous$M || frac != previous$frac
  rule = 'M from the held-out log-likelihood grid; frac = highest mean marker AUROC on a random half of the cells (seed 1), ties -> smaller frac; AUROC also reported on the other half. frac is tuned for marker agreement.'
  print(as.data.frame(summary))
  message(sprintf('%s: M = %s, frac = %s by marker AUROC (selection %.4f, held-out %.4f; %s)%s',
                  species, M, frac, summary$mean_auroc_selection[chosen_i],
                  summary$mean_auroc_heldout[chosen_i],
                  if (changed) sprintf('changed from M = %s, frac = %s', previous$M, previous$frac) else 'unchanged',
                  if (frac == max(fracs)) ' -- the largest frac won: the choice sits at the edge of the grid' else ''))

  out = file.path('data/02_annotated_data', paste0(species, '_ebcc_frac_auroc'))
  write_csv(summary, paste0(out, '.csv'))
  write_csv(per_class, paste0(out, '_per_class.csv'))
  write_csv(tibble(species = species, M = M, frac = frac, changed = changed,
                   previous_frac = previous$frac, rule = rule), paste0(out, '_chosen.csv'))

  d = summary |>
    select(frac, selection = mean_auroc_selection, heldout = mean_auroc_heldout) |>
    pivot_longer(-frac, names_to = 'half', values_to = 'mean_auroc')
  p = ggplot(d, aes(factor(frac), mean_auroc, colour = half, group = half)) +
    geom_line() + geom_point(size = 2.5) +
    geom_point(data = filter(d, frac == !!frac, half == 'selection'), shape = 21, size = 6,
               stroke = 1.4, colour = 'black') +
    scale_colour_manual(values = c(selection = '#0072B2', heldout = '#E69F00'),
                        labels = c(selection = 'Selection half', heldout = 'Held-out half'), name = NULL) +
    labs(title = sprintf('%s: EBCC frac by marker AUROC (M = %s)', species, M),
         subtitle = sprintf('circled: chosen frac = %s (frac tuned for marker agreement)', frac),
         x = 'frac (trusted-class prior strength)', y = 'Mean marker AUROC over called types') +
    theme_classic(base_size = 13)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(fig_dir, sprintf('ebcc_frac_auroc_%s.png', species)), p,
         width = 8, height = 5, dpi = 300, bg = 'white')
  list(frac = frac, changed = changed, fit = fits[[chosen_i]])
}
