# Section 2.2: are the differences between methods meaningful?  Rscript analysis/mcb_analysis.R
# Overall improvement over unintegrated (analysis/benchmark_figures.R output), dissected, one value per method x tissue.
# Main set (main_without_graph = TRUE): the 19 embedding methods x 4 tissues, graph-based methods (graph_only) excluded,
# for sections 1-8; all 23 methods as a sensitivity check (5a) and in mcb_results<suffix>_all23.csv.
# 1. lm(delta ~ method + tissue), tissue as a fixed block; residual = method x tissue. ANOVA, SS shares, partial eta^2.
# 2. Hsu's multiple comparisons with the best (MCB, Hsu 1996): D_i = mean_i - max_{j != i} mean_j,
#    interval [min(0, D_i - d SE), max(0, D_i + d SE)], SE = sqrt(2 MSE / n_tissue), d = one-sided Dunnett critical value
#    (k - 1 comparisons, residual df, correlation 0.5). 'Could be best' = upper > 0.
# 3. Practical equivalence: SD of paired differences between methods over the bootstrap subsamples (
#    all 19 embedding methods per tissue, as the additive model; graph-only methods have no valid bootstrap). Within delta of the best = MCB lower >= -delta.
# 4. Diagnostics (appendix). 5. Sensitivity (printout). 6. Main figure. 7. Friedman + Nemenyi critical difference (printout).
# Run for each score definition in `suffixes` ('' = reference, '_MAP_SAS' = with MAP / SAS), as analysis/benchmark_figures.R.
# Output (figures/mcb/): mcb_main<suffix>.{pdf,png}, mcb_diagnostics<suffix>.{pdf,png}, mcb_results<suffix>.csv; printout.
# Sections 2.3-2.5 (same outcome: improvement over unintegrated; tissue fixed; algorithms random so that conclusions are
# about strategies / inputs, not the particular algorithms; lmerTest, Kenward-Roger df, type II tests):
# 2.3 strategy: primary = native input only, y ~ tissue * strategy + (1 | algorithm) (residual = algorithm x tissue);
#     check = all 23 methods, + input as an additive fixed effect, + (1 | algorithm:tissue); also without input.
# 2.4 input: matched algorithms (incl. graph-only BBKNN; checks without it and without FastMNN) x 3 inputs x 4 tissues (balanced),
#     y ~ input * tissue + (1 | algorithm) + (1 | algorithm:input) for overall, bio and batch (type III); each input vs native
#     (emmeans) averaged over tissues and within each tissue. Figure: A raw gains (overall), B model by tissue, C-D raw gains
#     (bio, batch), E model by score; table grouping_input_model<suffix>.csv.
# 2.5 strategy x input: algorithms run on more than one input; grid of all 23 methods by strategy and input.
# Output (figures/mcb/): grouping_strategy<suffix>, grouping_input<suffix>, grouping_input_noBBKNN<suffix>, grouping_grid<suffix> .{pdf,png};
# grouping_grid<suffix>.csv, grouping_input_paired<suffix>.csv, grouping_input_model<suffix>.csv.
# needs only CRAN stats packages (no Seurat): run with R 4.5 as
#   Rscript --vanilla analysis/mcb_analysis.R
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr); library(stringr); library(ggplot2); library(patchwork)})

# ---------------- config ----------------
in_dir      = 'figures/benchmark_summary'
out         = 'figures/mcb'
suffixes    = c('_MAP_SAS', '')
alpha       = 0.05
deltas      = c(0.01, 0.02)                    # plus the median bootstrap SD of paired differences (added below)
delta_plot  = 0.01                              # placeholder band in the figure
top_n_boot  = 19                                # methods per tissue for the bootstrap SD of paired differences: 19 = all embedding methods
graph_only  = c('BBKNN', 'BBKNN UINMF', 'BBKNN scVI-latent', 'SAMap')
main_without_graph = TRUE                       # sections 1-8 on the 19 embedding methods; FALSE = all 23 (previous main)
label_based = 'scANVI'
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI')
input_levels = c('Native', 'scVI-latent', 'UINMF')
input_hl     = c(Native = '#d9d9d9', `scVI-latent` = '#7B5EA7', UINMF = '#CCBB44')   # as analysis/benchmark_figures.R
input_text   = c(Native = 'N', `scVI-latent` = 'S', UINMF = 'U')
ss_cols      = c(Method = '#3f007d', Tissue = '#E69F00', `Method x tissue` = '#bdbdbd')
# sections 2.3-2.5
strategy     = c(FastMNN = 'match', Scanorama = 'match', CCA = 'match', RPCA = 'match', SAMap = 'match',
                 BBKNN = 'balance', Harmony = 'balance', scVI = 'joint', scANVI = 'joint', LIGER = 'joint')   # as analysis/benchmark_figures.R
strategy_labels = c(match = 'Matching corresponding cells', balance = 'Balancing species in a shared space',
                    joint = 'Learning a joint representation')
strategy_cols   = c(match = '#4477AA', balance = '#228833', joint = '#AA3377')
matched_algorithms = c('Harmony', 'BBKNN', 'FastMNN', 'RPCA', 'CCA', 'Scanorama')   # run on all three inputs (2.4); BBKNN graph-only
embedding_matched  = setdiff(matched_algorithms, 'BBKNN')                              # 2.4 check: embedding methods only
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
tissue_cols   = c(bone_marrow = '#E69F00', colon = '#009E73', breast_cancer = '#CC79A7', fetal_liver = '#56B4E9')   # as analysis/benchmark_figures.R
tissue_shapes = c(bone_marrow = 16, colon = 17, breast_cancer = 15, fetal_liver = 18)
# -----------------------------------------
dir.create(out, showWarnings = FALSE, recursive = TRUE)
# Weighted Kendall tau (Vigna 2015), as analysis/benchmark_figures.R: hyperbolic weight 1/(r + 1) at importance rank r (0 = best),
# pair weight = sum of the two, averaged over both rankings as the importance order; = scipy.stats.weightedtau defaults.
wtau = function(x, y) {
  n = length(x)
  one = function(a, b) {
    rho = integer(n); rho[order(-a, -b)] = 0:(n - 1)
    w = 1 / (rho + 1); W = outer(w, w, '+'); up = upper.tri(W)
    sa = sign(outer(a, a, '-')); sb = sign(outer(b, b, '-'))
    sum((W * sa * sb)[up]) / sqrt(sum((W * (sa != 0))[up]) * sum((W * (sb != 0))[up]))
  }
  (one(x, y) + one(y, x)) / 2
}
label_of = function(m) {
  m = str_replace(m, '^(CCA|RPCA)', '\\1 (Seurat)')
  paste0(m, ifelse(m == label_based, ' †', ''), ifelse(str_remove(m, ' \\(Seurat\\)') %in% graph_only, ' ‡', ''))
}
input_of = function(m) factor(case_when(str_detect(m, 'UINMF') ~ 'UINMF', str_detect(m, 'scVI-latent') ~ 'scVI-latent',
                                        TRUE ~ 'Native'), levels = input_levels)

# ---- linear model + Hsu's MCB for one outcome ----
mcb = function(d, outcome) {
  d = d |> mutate(y = .data[[outcome]], method = factor(method), tissue = factor(tissue))
  stopifnot(nrow(d) == nlevels(d$method) * nlevels(d$tissue), all(table(d$method, d$tissue) == 1))
  fit = lm(y ~ method + tissue, data = d)
  a1 = anova(fit); a2 = anova(lm(y ~ tissue + method, data = d))
  k = nlevels(d$method); n_t = nlevels(d$tissue)
  df_res = fit$df.residual; mse = sum(resid(fit)^2) / df_res
  ss = c(Method = a1['method', 'Sum Sq'], Tissue = a1['tissue', 'Sum Sq'], `Method x tissue` = a1['Residuals', 'Sum Sq'])
  se = sqrt(2 * mse / n_t)
  R = matrix(0.5, k - 1, k - 1); diag(R) = 1
  set.seed(1)
  dcrit = mvtnorm::qmvt(1 - alpha, tail = 'lower.tail', df = df_res, corr = R)$quantile
  m = d |> group_by(method) |> summarise(mean = mean(y), .groups = 'drop') |> mutate(method = as.character(method))
  m$D = sapply(seq_len(nrow(m)), \(i) m$mean[i] - max(m$mean[-i]))
  tq = qt(1 - alpha / 2, df_res)
  m = m |> mutate(mcb_lo = pmin(0, D - dcrit * se), mcb_hi = pmax(0, D + dcrit * se),
                  ci_lo = D - tq * se, ci_hi = D + tq * se, could_be_best = mcb_hi > 0) |> arrange(-mean)
  list(fit = fit, anova = a1, anova_rev = a2, ss = ss, df_res = df_res, mse = mse, se = se, d = dcrit, res = m,
       partial_eta2 = ss[['Method']] / (ss[['Method']] + ss[['Method x tissue']]))
}
show_mcb = function(r, name) {
  g = r$res |> filter(could_be_best)
  cat(sprintf('%-34s d = %.3f, SE = %.4f, d*SE = %.4f; could be best: %d -> %s\n', name, r$d, r$se, r$d * r$se,
              nrow(g), paste(g$method, collapse = ', ')))
}

# ---- bootstrap: per-subsample scores (dissected), for paired differences between methods ----
boot_scores = function(extra_bio, extra_batch) {
  b = read_csv('bootstrap/summary/1_long.csv', show_col_types = FALSE) |> filter(`repeat` != 'point') |>
    transmute(tissue = dataset, rep = as.integer(`repeat`), method, metric, value = score)
  e = read_csv('bootstrap_extra/extra_long.csv', show_col_types = FALSE) |>
    pivot_longer(c(MAP, NC, SAS), names_to = 'metric', values_to = 'value') |> select(tissue, rep, method, metric, value)
  bind_rows(b, e) |> filter(metric %in% c(bio_metrics, batch_metrics, extra_bio, extra_batch),
                            !str_detect(method, '^Unintegrated'), !method %in% graph_only) |>
    mutate(group = ifelse(metric %in% c(batch_metrics, extra_batch), 'batch', 'bio')) |>
    group_by(tissue, rep, method, group) |> summarise(v = mean(value), n = n(), .groups = 'drop') |>
    filter(n == ifelse(group == 'bio', length(bio_metrics) + length(extra_bio), length(batch_metrics) + length(extra_batch))) |>
    select(-n) |> pivot_wider(names_from = group, values_from = v) |> mutate(overall = 0.6 * bio + 0.4 * batch) |>
    filter(!is.na(overall))
}

# ---- helpers for the figure of section 2.2 and sections 2.3-2.5 ----
alg_of      = function(m) word(m, 1)
strategy_of = function(m) factor(strategy[alg_of(m)], levels = names(strategy_labels))
fit_mm = function(f, data) suppressMessages(lmerTest::lmer(f, data = data))
kr = function(fit) as.data.frame(anova(fit, type = 2, ddf = 'Kenward-Roger'))
show_mm = function(fit, title) {
  cat(title, '\n')
  a = kr(fit)
  a |> transmute(term = rownames(a), NumDF, DenDF = round(DenDF, 1), F = round(`F value`, 2), p = signif(`Pr(>F)`, 3)) |>
    as.data.frame() |> print(row.names = FALSE)
  vc = as.data.frame(lme4::VarCorr(fit))
  cat('  random-effect SDs: ', paste(sprintf('%s %.4f', vc$grp, vc$sdcor), collapse = ', '),
      if (lme4::isSingular(fit)) ' (singular: a variance estimated at 0)' else '', '\n')
  invisible(a)
}
one_line = function(fit, terms) {   # compact F tests for the sensitivity printout
  a = kr(fit)
  paste(sapply(terms, \(t) sprintf('%s F(%d,%.1f) = %.2f, p = %.3f', t, a[t, 'NumDF'], a[t, 'DenDF'], a[t, 'F value'], a[t, 'Pr(>F)'])),
        collapse = '; ')
}
shares = function(fit, labels) { s = anova(fit)[['Sum Sq']]; setNames(s / sum(s), labels) }
show_shares = function(s) cat(paste(sprintf('  %-42s %5.1f%%', names(s), 100 * s), collapse = '\n'), '\n')
share_plot = function(s, highlight, xlab) {
  tibble(part = factor(names(s), levels = rev(names(s))), share = as.numeric(s), role = highlight[names(s)]) |>
    ggplot(aes(share, part, fill = role)) + geom_col(width = 0.7) +
    geom_text(aes(label = sprintf('%.0f%%', 100 * share)), hjust = -0.15, size = 3.2) +
    scale_fill_manual(values = c(group = '#3f007d', inter = '#b2182b', tissue = '#E69F00', other = '#bdbdbd'), guide = 'none') +
    scale_x_continuous(labels = scales::percent, expand = expansion(c(0, 0.15))) +
    labs(x = xlab, y = NULL) + theme_bw(base_size = 10.5) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.title.position = 'plot',
          plot.title = element_text(face = 'bold', size = 11))
}
vc_random = function(f, data, labels) {   # descriptive variance components: every term random, REML; SD per single value
  v = as.data.frame(lme4::VarCorr(suppressMessages(lme4::lmer(f, data = data, REML = TRUE))))
  setNames(v$sdcor[match(names(labels), v$grp)], labels)
}
sd_plot = function(sds, roles) {
  tibble(part = factor(names(sds), levels = rev(names(sds))), sd = as.numeric(sds), role = roles[names(sds)]) |>
    ggplot(aes(sd, part, fill = role)) + geom_col(width = 0.7) +
    geom_text(aes(label = sprintf('%.3f', sd)), hjust = -0.15, size = 3.2) +
    scale_fill_manual(values = c(group = '#3f007d', inter = '#b2182b', tissue = '#E69F00', other = '#bdbdbd'), guide = 'none') +
    scale_x_continuous(expand = expansion(c(0, 0.18))) +
    labs(x = 'Standard deviation (variance component)', y = NULL) + theme_bw(base_size = 10.5) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
}
paired = function(data, level, y = 'overall') data |> filter(input %in% c('Native', level)) |>
  select(algorithm, tissue, input, v = all_of(y)) |> pivot_wider(names_from = input, values_from = v) |>
  filter(!is.na(Native), !is.na(.data[[level]])) |> mutate(diff = .data[[level]] - Native, contrast = level)

for (suffix in suffixes) {
  extra_bio = if (suffix == '_MAP_SAS') 'MAP' else character(0); extra_batch = if (suffix == '_MAP_SAS') 'SAS' else character(0)
  cat(sprintf('\n\n################ Score definition: %s ################\n',
              if (suffix == '') 'reference (8 scib metrics)' else '+MAP (bio) +SAS (batch)'))
  d = read_csv(file.path(in_dir, paste0('benchmark_improvement', suffix, '.csv')), show_col_types = FALSE) |>
    select(method, tissue, overall = imp_overall, bio = imp_bio, batch = imp_batch)
  cat(sprintf('Design: %d rows, %d methods x %d tissues, every cell once: %s\n', nrow(d), n_distinct(d$method),
              n_distinct(d$tissue), all(table(d$method, d$tissue) == 1)))
  if (nrow(d) != 92 || !all(table(d$method, d$tissue) == 1)) stop('incomplete design: stop and report')
  d_all = d                                                   # all 23 methods (sensitivity 5a, mcb_results<suffix>_all23.csv)
  if (main_without_graph) d = filter(d_all, !method %in% graph_only)
  stopifnot(n_distinct(d$method) == if (main_without_graph) 19 else 23)
  main_name = sprintf('%d methods%s', n_distinct(d$method), if (main_without_graph) ', without graph-based' else '')

  # 1. model
  r = mcb(d, 'overall')
  cat(sprintf('\n== 1. ANOVA: lm(overall improvement ~ method + tissue), %s ==\n', main_name)); print(r$anova)
  cat(sprintf('Balanced check: SS(method) %.6g with method first, %.6g with tissue first (sequential = type III when equal)\n',
              r$anova['method', 'Sum Sq'], r$anova_rev['method', 'Sum Sq']))
  print(car::Anova(r$fit, type = 3))
  sh = r$ss / sum(r$ss)
  cat(sprintf('SS shares: method %.1f%%, tissue %.1f%%, method x tissue (residual) %.1f%%; partial eta^2 (method) = %.3f\n',
              100 * sh[1], 100 * sh[2], 100 * sh[3], r$partial_eta2))
  cat(sprintf('Residual df = %d, MSE = %.3g (residual SD %.4f)\n', r$df_res, r$mse, sqrt(r$mse)))

  # 2. MCB
  cat(sprintf('\n== 2. Hsu MCB, alpha = %.2f: d = %.4f (one-sided Dunnett, %d comparisons, df %d, corr 0.5), SE = %.4f ==\n',
              alpha, r$d, nrow(r$res) - 1, r$df_res, r$se))
  r$res |> transmute(method, mean = round(mean, 4), D = round(D, 4), mcb_lo = round(mcb_lo, 4), mcb_hi = round(mcb_hi, 4),
                     ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4), could_be_best) |> as.data.frame() |> print(row.names = FALSE)
  show_mcb(r, 'Could-be-best group (overall)')

  # 3. practical equivalence from the bootstrap
  bs = boot_scores(extra_bio, extra_batch)
  top = d |> filter(!method %in% graph_only) |> group_by(tissue) |> slice_max(overall, n = top_n_boot) |> ungroup()
  if (top_n_boot >= 19) stopifnot(all(count(top, tissue)$n == 19))   # every embedding method, every tissue
  pairs = top |> select(tissue, a = method) |> inner_join(top |> select(tissue, b = method), by = 'tissue',
                                                          relationship = 'many-to-many') |> filter(a < b)
  w = bs |> select(tissue, rep, method, overall)
  psd = pairs |> inner_join(w |> rename(a = method, ya = overall), by = c('tissue', 'a'), relationship = 'many-to-many') |>
    inner_join(w |> rename(b = method, yb = overall), by = c('tissue', 'b', 'rep')) |>
    group_by(tissue, a, b) |> summarise(sd = sd(ya - yb), n = n(), .groups = 'drop')
  delta_boot = median(psd$sd)
  cat(sprintf('\n== 3. Bootstrap SD of paired differences in overall improvement (top %d embedding methods per tissue, %d pairs, %d-%d subsamples) ==\n',
              top_n_boot, nrow(psd), min(psd$n), max(psd$n)))
  cat(sprintf('median %.4f, range %.4f-%.4f\n', delta_boot, min(psd$sd), max(psd$sd)))
  psd |> group_by(tissue) |> summarise(median = round(median(sd), 4), min = round(min(sd), 4), max = round(max(sd), 4)) |>
    as.data.frame() |> print(row.names = FALSE)
  all_d = c(deltas, delta_boot)
  names(all_d) = c(sprintf('within_%.2f', deltas), sprintf('within_bootSD_%.4f', delta_boot))
  for (nm in names(all_d)) r$res[[nm]] = r$res$could_be_best & r$res$mcb_lo >= -all_d[[nm]]
  cat('\nCould-be-best methods within delta of the best (MCB lower bound >= -delta):\n')
  r$res |> filter(could_be_best) |> select(method, mcb_lo, all_of(names(all_d))) |> mutate(mcb_lo = round(mcb_lo, 4)) |>
    as.data.frame() |> print(row.names = FALSE)

  # 4. diagnostics, on the model of the section 2.2 figure and ANOVA table (19 embedding methods, graph-based excluded)
  rdg = mcb(filter(d, !method %in% graph_only), 'overall')
  dg = tibble(method = as.character(rdg$fit$model$method), fitted = fitted(rdg$fit), resid = resid(rdg$fit))
  rsd = dg |> group_by(method) |> summarise(sd = sd(resid)) |> mutate(flag = sd > 2 * median(sd)) |> arrange(-sd)
  cat(sprintf('\n== 4. Diagnostics, %d methods (without graph-based). Residual SD by method: median %.4f; flagged (> 2x median): %s ==\n',
              n_distinct(dg$method), median(rsd$sd),
              if (any(rsd$flag)) paste(sprintf('%s (%.4f)', rsd$method[rsd$flag], rsd$sd[rsd$flag]), collapse = ', ') else 'none'))
  cat(sprintf('Shapiro-Wilk on residuals: W = %.3f, p = %.3f\n', shapiro.test(dg$resid)$statistic, shapiro.test(dg$resid)$p.value))
  th = theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
  q1 = ggplot(dg, aes(fitted, resid)) + geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey50') +
    geom_point(size = 1.6, alpha = 0.8) + geom_smooth(method = 'loess', formula = y ~ x, se = FALSE, colour = '#b2182b', linewidth = 0.6) +
    labs(x = 'Fitted', y = 'Residual', title = 'Residuals vs fitted') + th
  q2 = ggplot(dg, aes(sample = resid)) + stat_qq(size = 1.6, alpha = 0.8) + stat_qq_line(colour = '#b2182b') +
    labs(x = 'Theoretical quantile', y = 'Residual', title = 'Normal QQ') + th
  q3 = ggplot(rsd, aes(sd, factor(label_of(method), levels = rev(label_of(method))))) +
    geom_vline(xintercept = c(1, 2) * median(rsd$sd), linetype = c('solid', 'dashed'), colour = 'grey50') +
    geom_point(aes(colour = flag), size = 2) + scale_colour_manual(values = c(`FALSE` = 'grey20', `TRUE` = '#b2182b'), guide = 'none') +
    labs(x = 'Residual SD (line: median; dashed: 2x median)', y = NULL, title = 'Residual SD by method') + th
  pd = (q1 + q2) / q3 + patchwork::plot_layout(heights = c(1, 1.6))
  for (ext in c('pdf', 'png'))
    ggsave(file.path(out, paste0('mcb_diagnostics', suffix, '.', ext)), pd, width = 9, height = 10, dpi = 300, bg = 'white',
           device = if (ext == 'pdf') cairo_pdf else 'png')

  # 5. sensitivity
  cat('\n== 5. Sensitivity: could-be-best group ==\n')
  show_mcb(r, sprintf('(main) overall, %s', main_name))
  r_all = mcb(d_all, 'overall')
  if (main_without_graph) show_mcb(r_all, '(a) overall, all 23 methods (with graph-based)') else
    show_mcb(mcb(filter(d, !method %in% graph_only), 'overall'), '(a) overall, without graph-based')
  show_mcb(mcb(filter(d, method != label_based), 'overall'), '(b) overall, without scANVI')
  show_mcb(mcb(d, 'bio'), '(c) bio-conservation')
  show_mcb(mcb(d, 'batch'), '(c) batch correction')

  # 7. Friedman + Nemenyi on per-tissue ranks
  rk = d |> group_by(tissue) |> mutate(rk = rank(-overall)) |> ungroup()
  fr = friedman.test(rk ~ method | tissue, data = rk)
  k = n_distinct(d$method); N = n_distinct(d$tissue)
  cd = qtukey(1 - alpha, k, Inf) / sqrt(2) * sqrt(k * (k + 1) / (6 * N))
  cat(sprintf('\n== 7. Friedman chi^2 = %.2f, df = %d, p = %.4f; Nemenyi CD (alpha %.2f, %d methods, %d tissues) = %.2f ranks; mean-rank range %.2f-%.2f ==\n',
              fr$statistic, fr$parameter, fr$p.value, alpha, k, N, cd,
              min(tapply(rk$rk, rk$method, mean)), max(tapply(rk$rk, rk$method, mean))))

  # results table
  write_csv(r$res |> mutate(input = input_of(method), graph_only = method %in% graph_only, .after = method),
            file.path(out, paste0('mcb_results', suffix, '.csv')))
  if (main_without_graph)
    write_csv(r_all$res |> mutate(input = input_of(method), graph_only = method %in% graph_only, .after = method),
              file.path(out, paste0('mcb_results', suffix, '_all23.csv')))

  # 6. main figure
  ord = r$res$method   # = benchmark_improvement order (mean overall improvement)
  f = r$res |> mutate(y = factor(label_of(method), levels = rev(label_of(ord))), input = input_of(method),
                      pos = seq_along(method))
  strip = ggplot(f, aes(1, y)) +
    geom_text(aes(x = 0, label = pos), size = 3.4, fontface = 'bold', colour = 'grey10') +
    geom_tile(aes(fill = input), width = 0.9, height = 0.85) +
    geom_text(aes(label = input_text[input], colour = input == 'Native'), size = 2.6, fontface = 'bold') +
    scale_x_continuous(limits = c(-0.5, 1.5)) +
    scale_fill_manual(values = input_hl, name = 'Input', drop = FALSE, labels = \(x) sprintf('%s (%s)', x, input_text[x])) +
    scale_colour_manual(values = c(`TRUE` = 'grey20', `FALSE` = 'white'), guide = 'none') +
    theme_void(base_size = 11) +
    theme(axis.text.y = element_text(colour = 'grey10', hjust = 1, size = 9.5, margin = margin(r = 2)))
  pa = ggplot(f, aes(D, y)) +
    geom_blank() +
    annotate('rect', xmin = -delta_plot, xmax = 0, ymin = -Inf, ymax = Inf, fill = '#7B5EA7', alpha = 0.12) +
    geom_vline(xintercept = 0, colour = 'grey30', linewidth = 0.5) +
    geom_errorbar(aes(xmin = mcb_lo, xmax = mcb_hi, colour = could_be_best), width = 0, orientation = 'y', linewidth = 0.9) +
    geom_point(aes(shape = could_be_best, colour = could_be_best), size = 2.6, fill = 'black', stroke = 0.9) +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), labels = c(`TRUE` = 'Could be best', `FALSE` = 'Worse than the best'),
                       name = NULL, breaks = c('TRUE', 'FALSE')) +
    scale_colour_manual(values = c(`TRUE` = 'black', `FALSE` = 'grey55'), labels = c(`TRUE` = 'Could be best',
                        `FALSE` = 'Worse than the best'), name = NULL, breaks = c('TRUE', 'FALSE')) +
    labs(x = 'Mean overall improvement minus the best of the other methods (Hsu MCB, 95%)', y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(colour = 'grey94'),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  ssd = tibble(part = names(sh), share = as.numeric(sh)) |> mutate(hi = cumsum(share), lo = hi - share)
  pb = ggplot(ssd) +
    geom_rect(aes(xmin = lo, xmax = hi, ymin = 0, ymax = 1, fill = part), colour = 'white', linewidth = 0.8) +
    geom_text(aes(x = (lo + hi) / 2, y = 0.5, label = sprintf('%s\n%.0f%%', part, 100 * share), colour = part == 'Method'),
              size = 3.2, fontface = 'bold', lineheight = 0.9) +
    scale_fill_manual(values = ss_cols, guide = 'none') +
    scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = 'grey10'), guide = 'none') +
    scale_x_continuous(labels = scales::percent, expand = expansion(0), limits = c(0, 1)) +
    labs(x = 'Share of the sum of squares of overall improvement', y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank())
  pm = (strip + pa + patchwork::plot_layout(widths = c(0.16, 1))) / pb +
    patchwork::plot_layout(heights = c(1, 0.12), guides = 'collect') +
    patchwork::plot_annotation(tag_levels = list(c('a', '', 'b'))) &
    theme(legend.position = 'bottom', plot.tag = element_text(face = 'bold', size = 14))
  for (ext in c('pdf', 'png'))
    ggsave(file.path(out, paste0('mcb_main', suffix, '.', ext)), pm, width = 10, height = 8.5, dpi = 300, bg = 'white',
           device = if (ext == 'pdf') cairo_pdf else 'png')

  # 8. figure for section 2.2: MCB (A), agreement between the tissues' rankings (B), variance shares (C); notes go in the caption
  # graph-based methods (graph_only) excluded: the model, MCB, tau and variance components are refitted on the embedding methods
  d22 = filter(d, !method %in% graph_only)
  r22 = mcb(d22, 'overall')
  cat(sprintf('\n== 8. Section 2.2 figure: %d methods (without graph-based) x %d tissues ==\n', n_distinct(d22$method), n_distinct(d22$tissue)))
  show_mcb(r22, 'Could-be-best group (overall)')
  rmat = d22 |> group_by(tissue) |> mutate(r = rank(-overall)) |> ungroup() |> select(method, tissue, r) |>
    pivot_wider(names_from = tissue, values_from = r)
  rmat = as.matrix(rmat[, names(tissue_labels)])
  # panel B: weighted Kendall tau (wtau, as analysis/benchmark_figures.R: Vigna 2015, hyperbolic weights = scipy.stats.weightedtau);
  # the plain Kendall tau is printed next to it
  tau = expand_grid(a = names(tissue_labels), b = names(tissue_labels)) |>
    mutate(tau = map2_dbl(a, b, \(x, y) cor(rmat[, x], rmat[, y], method = 'kendall')),
           wtau = map2_dbl(a, b, \(x, y) wtau(-rmat[, x], -rmat[, y])),   # rank 1 = best = most important
           ia = match(a, names(tissue_labels)), ib = match(b, names(tissue_labels)))
  cat('Kendall tau and weighted Kendall tau between tissues (without graph-based):\n')
  tau |> filter(ia < ib) |> transmute(pair = paste(a, b, sep = ' vs '), tau = round(tau, 2), weighted_tau = round(wtau, 2)) |>
    as.data.frame() |> print(row.names = FALSE)
  write_csv(tau |> filter(ia != ib) |> select(a, b, tau, weighted_tau = wtau), file.path(out, paste0('section22_tissue_tau', suffix, '.csv')))
  # Kendall's W over the four tissues on the same 19 methods, permutation p as analysis/benchmark_figures.R (ranks shuffled within tissue)
  kendall_w = function(m) { k = ncol(m); n = nrow(m); R = rowSums(m); 12 * sum((R - mean(R))^2) / (k^2 * (n^3 - n)) }
  W22 = kendall_w(rmat); set.seed(1)
  W22_p = (1 + sum(replicate(10000, kendall_w(apply(rmat, 2, sample))) >= W22)) / (10000 + 1)
  cat(sprintf("Kendall's W (without graph-based, %d methods) = %.2f, permutation p = %.4f (10,000 permutations)\n", nrow(rmat), W22, W22_p))
  write_csv(tibble(n_methods = nrow(rmat), W = W22, perm_p = W22_p), file.path(out, paste0('section22_kendall_w', suffix, '.csv')))
  tl2 = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
  tau_off = filter(tau, ia > ib)
  pt2 = ggplot(tau, aes(factor(tl2[b], levels = tl2), factor(tl2[a], levels = rev(tl2)))) +
    geom_tile(data = filter(tau, ia == ib), fill = 'grey85', colour = 'white', linewidth = 1) +
    geom_text(data = filter(tau, ia == ib), label = '1', colour = 'grey45', size = 3.8) +
    geom_tile(data = tau_off, aes(fill = wtau), colour = 'white', linewidth = 1) +
    geom_text(data = tau_off, aes(label = sprintf('%.2f', wtau), colour = abs(wtau) > 0.5), size = 3.8, fontface = 'bold') +
    scale_fill_gradient2(low = '#b2182b', mid = '#f7f7f7', high = '#3f007d', midpoint = 0, limits = c(-1, 1), name = "Weighted Kendall's τ",
                         guide = guide_colourbar(title.vjust = 0.8, barwidth = 7, barheight = 0.6)) +
    scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = 'grey15'), guide = 'none') +
    scale_x_discrete(labels = \(x) str_replace(x, ' ', '\n')) +
    coord_equal() + labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10.5) + theme(panel.grid = element_blank(), axis.text = element_text(colour = 'grey10'))
  # variance components (ANOVA estimators; per single method x tissue value), and cell sampling on the same scale:
  # SD of a paired difference between two methods over the bootstrap / sqrt(2)
  a = r22$anova; n_m = nlevels(factor(d22$method)); n_t = nlevels(factor(d22$tissue))
  vc = c(Method = (a['method', 'Mean Sq'] - r22$mse) / n_t, Tissue = (a['tissue', 'Mean Sq'] - r22$mse) / n_m,
         `Method x tissue` = r22$mse, `Cell sampling` = (delta_boot / sqrt(2))^2)
  vc_sd = sqrt(pmax(vc, 0))
  cat(sprintf('\n== 8b. Variance components (SD): method %.4f, tissue %.4f, method x tissue %.4f, cell sampling %.4f ==\n',
              vc_sd[1], vc_sd[2], vc_sd[3], vc_sd[4]))
  cat(sprintf('method x tissue / method = %.2f; method x tissue / cell sampling = %.1f (same ratio on the paired-difference scale: %.4f vs %.4f)\n',
              vc_sd[3] / vc_sd[1], vc_sd[3] / vc_sd[4], sqrt(2) * vc_sd[3], delta_boot))
  pv2 = tibble(part = factor(names(vc_sd), levels = rev(names(vc_sd))), sd = vc_sd,
               role = c('group', 'tissue', 'inter', 'other')) |>
    ggplot(aes(sd, part, fill = role)) + geom_col(width = 0.7) +
    geom_text(aes(label = sprintf('%.3f', sd)), hjust = -0.15, size = 3.2) +
    scale_fill_manual(values = c(group = '#3f007d', inter = '#b2182b', tissue = '#E69F00', other = '#bdbdbd'), guide = 'none') +
    scale_x_continuous(expand = expansion(c(0, 0.18))) +
    labs(x = 'Standard deviation (variance component)', y = NULL) + theme_bw(base_size = 10.5) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
  f22 = r22$res |> mutate(y = factor(label_of(method), levels = rev(label_of(r22$res$method))), input = input_of(method), pos = seq_along(method))
  strip22 = strip %+% f22
  pa2 = pa %+% f22; pa2$layers = pa2$layers[-2]   # without the placeholder equivalence band
  pa2 = pa2 + labs(x = 'Mean overall improvement minus the best of the other methods\n(Hsu MCB, 95% simultaneous)')
  p22 = (strip22 + pa2 + patchwork::plot_layout(widths = c(0.16, 1))) - (pt2 / pv2 + patchwork::plot_layout(heights = c(1, 0.5))) +
    patchwork::plot_layout(widths = c(1.45, 1), guides = 'collect') +
    patchwork::plot_annotation(tag_levels = list(c('A', '', 'B', 'C'))) &
    theme(legend.position = 'bottom', plot.tag = element_text(face = 'bold', size = 14))
  for (ext in c('pdf', 'png'))
    ggsave(file.path(out, paste0('section22', suffix, '.', ext)), p22, width = 13, height = 8.5, dpi = 300, bg = 'white',
           device = if (ext == 'pdf') cairo_pdf else 'png')
}

for (suffix in suffixes) {
  cat(sprintf('\n\n################ Sections 2.3-2.5, score definition: %s ################\n',
              if (suffix == '') 'reference (8 scib metrics)' else '+MAP (bio) +SAS (batch)'))
  op = options(contrasts = c('contr.sum', 'contr.poly'))
  d = read_csv(file.path(in_dir, paste0('benchmark_improvement', suffix, '.csv')), show_col_types = FALSE) |>
    transmute(method, tissue = factor(tissue, levels = names(tissue_labels)), overall = imp_overall, bio = imp_bio, batch = imp_batch,
              algorithm = alg_of(method), input = input_of(method), strategy = strategy_of(method))
  stopifnot(!anyNA(d$strategy), nrow(d) == 92)

  # ---- 2.3 strategy ----
  nat = filter(d, input == 'Native')
  stopifnot(all(table(nat$algorithm, nat$tissue) == 1))
  cat(sprintf('\n== 2.3 Strategy. Native input: %d algorithms (%s) x %d tissues ==\n', n_distinct(nat$algorithm),
              paste(sprintf('%s %d', names(strategy_labels), as.integer(table(strategy_of(unique(nat$algorithm))))), collapse = ', '),
              nlevels(nat$tissue)))
  s23 = shares(lm(overall ~ tissue + strategy + strategy:algorithm + strategy:tissue, data = nat),
               c('Tissue', 'Strategy', 'Algorithm within strategy', 'Strategy x tissue', 'Algorithm x tissue (residual)'))
  cat('SS shares (native input, nested, orthogonal):\n'); show_shares(s23)
  ms = anova(lm(overall ~ tissue + strategy + strategy:algorithm + strategy:tissue, data = nat))
  cat(sprintf('  mean squares: strategy %.5f, algorithm within strategy %.5f, strategy x tissue %.5f, algorithm x tissue %.5f\n',
              ms[2, 'Mean Sq'], ms[3, 'Mean Sq'], ms[4, 'Mean Sq'], ms[5, 'Mean Sq']))
  v23 = vc_random(overall ~ 1 + (1 | tissue) + (1 | strategy) + (1 | algorithm) + (1 | strategy:tissue), nat,
                  c(strategy = 'Strategy', algorithm = 'Algorithm within strategy', tissue = 'Tissue', `strategy:tissue` = 'Strategy x tissue',
                    Residual = 'Algorithm x tissue (residual)'))
  cat('Variance components, SD (native input, all terms random, REML):\n'); cat(sprintf('  %-32s %.4f\n', names(v23), v23), sep = '')
  f_nat = fit_mm(overall ~ tissue * strategy + (1 | algorithm), nat)
  show_mm(f_nat, '(primary) native input: overall ~ tissue * strategy + (1 | algorithm)')
  f_all = fit_mm(overall ~ tissue * strategy + (1 | algorithm) + (1 | algorithm:tissue), d)
  show_mm(f_all, '(check) all 23 methods, input ignored: overall ~ tissue * strategy + (1 | algorithm) + (1 | algorithm:tissue)')
  f_adj = fit_mm(overall ~ tissue * strategy + input + (1 | algorithm) + (1 | algorithm:tissue), d)
  show_mm(f_adj, '(check) all 23 methods, input additive: overall ~ tissue * strategy + input + (1 | algorithm) + (1 | algorithm:tissue)')
  em = bind_rows(lapply(list(`Native input only` = f_nat, `All methods, input ignored` = f_all,
                             `All methods, adjusted for input` = f_adj), \(f)
    suppressMessages(as.data.frame(emmeans::emmeans(f, ~ strategy, lmer.df = 'kenward-roger')))), .id = 'model') |>
    mutate(model = factor(model, levels = c('Native input only', 'All methods, input ignored', 'All methods, adjusted for input')))
  cat('Strategy means (emmeans, averaged over tissues), 95% CI:\n')
  em |> transmute(model, strategy, mean = round(emmean, 4), lo = round(lower.CL, 4), hi = round(upper.CL, 4), df = round(df, 1)) |>
    as.data.frame() |> print(row.names = FALSE)
  pt = nat |> group_by(tissue, strategy) |> summarise(m = mean(overall), .groups = 'drop')
  cat('Native input, strategy mean per tissue (best marked *):\n')
  pt |> group_by(tissue) |> mutate(txt = sprintf('%.3f%s', m, ifelse(m == max(m), '*', ''))) |> select(-m) |>
    pivot_wider(names_from = strategy, values_from = txt) |> as.data.frame() |> print(row.names = FALSE)
  cat(sprintf('Algorithm fixed instead (strategy tested against algorithm x tissue): strategy F(%d,%d) = %.2f, p = %.3f; strategy x tissue F(%d,%d) = %.2f, p = %.3f\n',
              ms[2, 'Df'], ms[5, 'Df'], ms[2, 'F value'], ms[2, 'Pr(>F)'], ms[4, 'Df'], ms[5, 'Df'], ms[4, 'F value'], ms[4, 'Pr(>F)']))
  cat('Sensitivity (native input, primary model):\n')
  ng = filter(nat, !method %in% graph_only)
  cat('  without graph-only (balance = Harmony only): ', one_line(fit_mm(overall ~ tissue * strategy + (1 | algorithm), ng), c('strategy', 'tissue:strategy')), '\n')
  for (y in c('bio', 'batch'))
    cat(sprintf('  %-5s: %s\n', y, one_line(fit_mm(reformulate(c('tissue * strategy', '(1 | algorithm)'), y), nat), c('strategy', 'tissue:strategy'))))
  cat('  all 23, input additive, without graph-only: ',
      one_line(fit_mm(overall ~ tissue * strategy + input + (1 | algorithm) + (1 | algorithm:tissue), filter(d, !method %in% graph_only)),
               c('strategy', 'tissue:strategy', 'input')), '\n')

  # ---- 2.4 input ----
  inp = filter(d, algorithm %in% matched_algorithms)
  stopifnot(nrow(inp) == length(matched_algorithms) * 3 * 4, all(table(inp$algorithm, inp$input, inp$tissue) == 1))
  cat(sprintf('\n== 2.4 Input. %d algorithms (%s) x 3 inputs x 4 tissues, balanced ==\n', length(matched_algorithms),
              paste(matched_algorithms, collapse = ', ')))
  s24 = shares(lm(overall ~ tissue + algorithm + input + tissue:algorithm + tissue:input + algorithm:input, data = inp),
               c('Tissue', 'Algorithm', 'Input', 'Algorithm x tissue', 'Input x tissue', 'Algorithm x input',
                 'Algorithm x input x tissue (residual)'))
  cat('SS shares (orthogonal):\n'); show_shares(s24)
  # model, fitted to the overall, bio and batch improvement: y ~ input * tissue + (1 | algorithm) + (1 | algorithm:input).
  # Input and tissue fixed (sum-to-zero); algorithms random, algorithm:input = how much an input's benefit varies between
  # algorithms (the error term for input). Type III F tests, Kenward-Roger df. Each input vs native (emmeans, trt.vs.ctrl),
  # averaged over tissues and within each tissue. Sets: all matched algorithms (main), without BBKNN, without FastMNN (checks).
  fit24 = function(data, y) suppressMessages(lmerTest::lmer(as.formula(paste(y, '~ input * tissue + (1 | algorithm) + (1 | algorithm:input)')), data = data))
  sets24 = list(`All matched algorithms` = matched_algorithms, `Without BBKNN (embedding only)` = embedding_matched,
                `Without FastMNN` = setdiff(matched_algorithms, 'FastMNN'))
  vs_native = function(em) as.data.frame(summary(emmeans::contrast(em, 'trt.vs.ctrl', adjust = 'none'), infer = TRUE))
  a24 = list(); c24 = list()
  for (st in names(sets24)) for (y in c('overall', 'bio', 'batch')) {
    m = fit24(filter(inp, algorithm %in% sets24[[st]]), y)
    a = anova(m, type = 3, ddf = 'Kenward-Roger'); vc = as.data.frame(lme4::VarCorr(m))
    a24[[paste(st, y)]] = as.data.frame(a) |> tibble::rownames_to_column('term') |>
      transmute(set = st, score = y, term, NumDF, DenDF, F = `F value`, p = `Pr(>F)`,
                random_sd = paste(sprintf('%s %.4f', vc$grp, vc$sdcor), collapse = ', '), singular = lme4::isSingular(m))
    c24[[paste(st, y)]] = bind_rows(vs_native(emmeans::emmeans(m, ~ input, lmer.df = 'kenward-roger')) |> mutate(tissue = 'all'),
                                    vs_native(emmeans::emmeans(m, ~ input | tissue, lmer.df = 'kenward-roger')) |> mutate(tissue = as.character(tissue))) |>
      transmute(set = st, score = y, tissue, contrast = str_remove(gsub('[()]', '', contrast), ' - Native'),
                est = estimate, lo = lower.CL, hi = upper.CL, df, p = p.value)
  }
  a24 = bind_rows(a24); c24 = bind_rows(c24)
  write_csv(bind_rows(a24 |> mutate(table = 'anova'), c24 |> mutate(table = 'vs_native')) |> relocate(table),
            file.path(out, paste0('grouping_input_model', suffix, '.csv')))
  cat('Model: y ~ input * tissue + (1 | algorithm) + (1 | algorithm:input); type III, Kenward-Roger\n')
  for (st in names(sets24)) for (y in c('overall', 'bio', 'batch')) {
    cat(sprintf('-- %s, %s --\n', st, y))
    a = filter(a24, set == st, score == y)
    a |> transmute(term, NumDF, DenDF = round(DenDF, 1), F = round(F, 2), p = signif(p, 3)) |> as.data.frame() |> print(row.names = FALSE)
    cat('  random-effect SDs: ', a$random_sd[1], if (a$singular[1]) ' (singular: a variance estimated at 0)' else '', '\n')
    cat('  each input vs native (95% CI):\n')
    filter(c24, set == st, score == y) |> transmute(tissue = ifelse(tissue == 'all', 'mean over tissues', tissue_labels[tissue]), contrast,
                                                    est = sprintf('%+.4f [%+.4f, %+.4f]', est, lo, hi), df = round(df, 1), p = signif(p, 3)) |>
      as.data.frame() |> print(row.names = FALSE)
  }
  pd = bind_rows(paired(inp, 'scVI-latent'), paired(inp, 'UINMF'))
  write_csv(pd, file.path(out, paste0('grouping_input_paired', suffix, '.csv')))
  # raw paired gains in the two components (figure panels D, E; descriptive)
  pdc = bind_rows(lapply(c('bio', 'batch'), \(y) bind_rows(paired(inp, 'scVI-latent', y), paired(inp, 'UINMF', y)) |> mutate(score = y)))

  # ---- 2.5 strategy x input ----
  multi = d |> group_by(algorithm) |> filter(n_distinct(input) > 1) |> ungroup()
  cat(sprintf('\n== 2.5 Strategy x input. Algorithms on more than one input: %s ==\n',
              paste(sprintf('%s (%s)', unique(multi$algorithm), strategy[unique(multi$algorithm)]), collapse = ', ')))
  show_mm(fit_mm(overall ~ tissue + strategy * input + tissue:strategy + tissue:input + (1 | algorithm) + (1 | algorithm:tissue) +
                   (1 | algorithm:input), multi),
          'overall ~ tissue + strategy * input + tissue:strategy + tissue:input + (1 | algorithm) + (1 | algorithm:tissue) + (1 | algorithm:input)')
  gain = bind_rows(paired(d, 'scVI-latent'), paired(d, 'UINMF')) |> mutate(strategy = strategy_of(algorithm), input = contrast)
  grid = d |> group_by(strategy, input) |>
    summarise(algorithms = paste(sort(unique(algorithm)), collapse = ', '), n_alg = n_distinct(algorithm), mean = mean(overall),
              worst_tissue = min(tapply(overall, tissue, mean)), tissues_pos = sum(tapply(overall, tissue, mean) > 0), .groups = 'drop') |>
    left_join(gain |> group_by(strategy, input, tissue) |> summarise(g = mean(diff), .groups = 'drop') |> group_by(strategy, input) |>
                summarise(gain_vs_native = mean(g), tissues_beat_native = sum(g > 0), .groups = 'drop'), by = c('strategy', 'input')) |>
    complete(strategy, input)
  write_csv(grid, file.path(out, paste0('grouping_grid', suffix, '.csv')))
  cat('Grid (all 23 methods): mean improvement over algorithms and tissues; worst tissue; gain over the same algorithms on native input:\n')
  grid |> mutate(across(where(is.double), \(x) round(x, 3))) |> as.data.frame() |> print(row.names = FALSE)
  options(op)

  # ---- figures ----
  th = theme_bw(base_size = 10.5) + theme(panel.grid.minor = element_blank(), plot.title.position = 'plot',
                                          plot.title = element_text(face = 'bold', size = 11))
  dodge = position_dodge(width = 0.55)
  pa = ggplot(nat, aes(tissue, overall, colour = strategy, group = strategy)) +
    geom_hline(yintercept = 0, colour = 'grey50', linetype = 'dashed') +
    geom_point(position = dodge, size = 1.7, alpha = 0.45) +
    geom_line(data = pt, aes(y = m), position = dodge, linewidth = 0.9) +
    geom_point(data = pt, aes(y = m, fill = strategy), position = dodge, shape = 23, size = 3.3, colour = 'black') +
    scale_x_discrete(labels = tissue_labels) +
    scale_colour_manual(values = strategy_cols, labels = strategy_labels, name = NULL) +
    scale_fill_manual(values = strategy_cols, labels = strategy_labels, name = NULL) +
    labs(x = NULL, y = 'Overall improvement over unintegrated') +
    th + theme(legend.position = 'bottom', legend.direction = 'vertical')
  pb = sd_plot(v23, c(Strategy = 'group', `Algorithm within strategy` = 'other', Tissue = 'tissue', `Strategy x tissue` = 'inter',
                      `Algorithm x tissue (residual)` = 'other'))
  pc = ggplot(em, aes(emmean, factor(strategy, levels = rev(names(strategy_labels)), labels = rev(strategy_labels)), colour = model)) +
    geom_vline(xintercept = 0, colour = 'grey50', linetype = 'dashed') +
    geom_pointrange(aes(xmin = lower.CL, xmax = upper.CL), position = position_dodge(width = -0.6), size = 0.4) +
    scale_colour_manual(values = c('#000000', '#999999', '#E69F00'), name = NULL) +
    labs(x = 'Strategy mean improvement (95% CI, algorithms random)', y = NULL) +
    th + theme(legend.position = 'bottom', legend.direction = 'vertical')
  p23 = (pa | (pb / pc + patchwork::plot_layout(heights = c(0.8, 1)))) + patchwork::plot_layout(widths = c(1.25, 1)) +
    patchwork::plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 13))
  # 2.4 figure: (A) raw paired gains (overall) by algorithm and tissue; (B) model: each input vs native within each tissue and
  # averaged over tissues (overall); (D, E) raw paired gains in bio-conservation and batch correction; (F) model, averaged over
  # tissues, by score (panels tagged A-E). For one set of matched algorithms.
  lab_c = c(`scVI-latent` = 'scVI-latent minus native', UINMF = 'UINMF minus native')
  input_fig = function(st, algs) {
  pq = pd |> filter(algorithm %in% algs) |> mutate(algorithm = factor(label_of(algorithm), levels = rev(label_of(algs))),
                    contrast = factor(lab_c[contrast], levels = lab_c)) |>
    ggplot(aes(diff, algorithm)) +
    geom_vline(xintercept = 0, colour = 'grey40') +
    geom_point(aes(colour = tissue, shape = tissue), size = 2.4, alpha = 0.9) +
    facet_wrap(~ contrast, nrow = 1) +
    scale_colour_manual(values = tissue_cols, labels = tissue_labels, name = NULL) +
    scale_shape_manual(values = tissue_shapes, labels = tissue_labels, name = NULL) +
    labs(x = 'Gain in overall improvement over the same algorithm on native input', y = NULL) +
    th + theme(legend.position = 'bottom', strip.background = element_rect(fill = 'grey95', colour = 'grey70'),
               strip.text = element_text(face = 'bold'))
  rows = c(tissue_labels, all = 'Mean over tissues')
  pB = c24 |> filter(set == st, score == 'overall') |>
    mutate(row = factor(rows[tissue], levels = rev(rows)), contrast = factor(lab_c[contrast], levels = lab_c)) |>
    ggplot(aes(est, row)) +
    geom_vline(xintercept = 0, colour = 'grey40') +
    geom_pointrange(aes(xmin = lo, xmax = hi, colour = tissue, shape = tissue), size = 0.5) +
    facet_wrap(~ contrast, ncol = 1) +
    scale_colour_manual(values = c(tissue_cols, all = 'black'), guide = 'none') +
    scale_shape_manual(values = c(tissue_shapes, all = 18), guide = 'none') +
    labs(x = 'Gain in overall improvement (model, 95% CI)', y = NULL) + theme_bw(base_size = 10.5) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = 'grey95', colour = 'grey70'), strip.text = element_text(face = 'bold'))
  comp_plot = function(y, xlab) pdc |> filter(score == y, algorithm %in% algs) |>
    mutate(algorithm = factor(label_of(algorithm), levels = rev(label_of(algs))), contrast = factor(lab_c[contrast], levels = lab_c)) |>
    ggplot(aes(diff, algorithm)) +
    geom_vline(xintercept = 0, colour = 'grey40') +
    geom_point(aes(colour = tissue, shape = tissue), size = 2.4, alpha = 0.9) +
    facet_wrap(~ contrast, nrow = 1) +
    scale_colour_manual(values = tissue_cols, labels = tissue_labels, guide = 'none') +
    scale_shape_manual(values = tissue_shapes, labels = tissue_labels, guide = 'none') +
    labs(x = xlab, y = NULL) +
    th + theme(strip.background = element_rect(fill = 'grey95', colour = 'grey70'), strip.text = element_text(face = 'bold'))
  pD = comp_plot('bio', 'Gain in bio-conservation improvement')
  pE = comp_plot('batch', 'Gain in batch-correction improvement')
  pF = c24 |> filter(set == st, tissue == 'all') |>
    mutate(score = factor(score, levels = c('batch', 'bio', 'overall'), labels = c('Batch correction', 'Bio-conservation', 'Overall')),
           contrast = factor(lab_c[contrast], levels = lab_c)) |>
    ggplot(aes(est, score)) +
    geom_vline(xintercept = 0, colour = 'grey40') +
    geom_pointrange(aes(xmin = lo, xmax = hi), shape = 18, size = 0.6) +
    facet_wrap(~ contrast, ncol = 1) +
    labs(x = 'Gain, mean over tissues (model, 95% CI)', y = NULL) + theme_bw(base_size = 10.5) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = 'grey95', colour = 'grey70'), strip.text = element_text(face = 'bold'))
  top = (pq | pB) + patchwork::plot_layout(widths = c(1.6, 1))
  bottom = (pD | pE | pF) + patchwork::plot_layout(widths = c(1, 1, 0.75))
  (top / bottom) + patchwork::plot_layout(heights = c(1, 0.8)) + patchwork::plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(face = 'bold', size = 13))
  }
  p24 = input_fig('All matched algorithms', matched_algorithms)
  p24e = input_fig('Without BBKNN (embedding only)', embedding_matched)   # check without the graph-only BBKNN
  gg = grid |> mutate(s = factor(strategy, levels = rev(names(strategy_labels)), labels = rev(strategy_labels)),
                      input = factor(input, levels = input_levels),
                      txt = ifelse(is.na(mean), 'not run', paste0(sprintf('%.3f', mean), '\nworst tissue ', sprintf('%.3f', worst_tissue),
                                   ifelse(is.na(tissues_beat_native), '', sprintf('\nbeats native in %d/4', tissues_beat_native)),
                                   '\n', str_replace_all(algorithms, '(([^,]+, ){2})', '\\1\n'))))
  pg = ggplot(gg, aes(input, s)) +
    geom_tile(aes(fill = mean), colour = 'white', linewidth = 1) +
    geom_text(aes(label = txt, colour = !is.na(mean) & mean > 0.14), size = 2.9, lineheight = 0.95) +
    scale_fill_gradient(low = '#f2f0f7', high = '#3f007d', na.value = 'grey92', name = 'Mean overall\nimprovement') +
    scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = 'grey10'), guide = 'none') +
    scale_x_discrete(position = 'top') + labs(x = 'Input', y = NULL) +
    theme_minimal(base_size = 11) + theme(panel.grid = element_blank(), axis.text = element_text(colour = 'grey10'))
  for (ext in c('pdf', 'png')) {
    dev = if (ext == 'pdf') cairo_pdf else 'png'
    ggsave(file.path(out, paste0('grouping_strategy', suffix, '.', ext)), p23, width = 13, height = 7, dpi = 300, bg = 'white', device = dev)
    ggsave(file.path(out, paste0('grouping_input', suffix, '.', ext)), p24, width = 15, height = 11.5, dpi = 300, bg = 'white', device = dev)
    ggsave(file.path(out, paste0('grouping_input_noBBKNN', suffix, '.', ext)), p24e, width = 15, height = 11.5, dpi = 300, bg = 'white', device = dev)
    ggsave(file.path(out, paste0('grouping_grid', suffix, '.', ext)), pg, width = 9, height = 5, dpi = 300, bg = 'white', device = dev)
  }
}
