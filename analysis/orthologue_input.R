# Does the input advantage hold as shared orthologues are removed?  Rscript analysis/orthologue_input.R
# Cases: dissected (0%) and genes2/4/6/8 (20/40/60/80% of shared one-to-one orthologues removed completely), 4 tissues.
# Scores exactly as analysis/benchmark_figures.R: bio = mean(NMI, ARI, ASW_label, cLISI [+ MAP]), batch = mean(PCR_batch,
# ASW_label/batch, iLISI, graph_conn [+ SAS]), overall = 0.6 bio + 0.4 batch.
# Gain = overall of an algorithm on scVI-latent (or UINMF) input minus the same algorithm on native input, in the same
# tissue and case (= the difference in improvement over unintegrated, since the reference cancels).
# Model, per input, on the matched algorithms (as section 2.4):
#   gain ~ algorithm * removed + (1 + removed || tissue) + (1 | tissue:case)
#   algorithm fixed (claim about these algorithms); removed = per 10% of shared orthologues removed (slope);
#   tissue random intercept and slope (generalise to new datasets); tissue:case = dataset, shared by the algorithms.
#   Sum-to-zero contrasts: intercept / slope = mean over the algorithms. Kenward-Roger df.
# Run for the reference definition and, once MAP / SAS exist for every reduced-gene case, for +MAP +SAS.
# Section 3.2.4 (Is UINMF the answer to low shared genes?), same model on paired differences:
#   UINMF vs native, UINMF vs scVI-latent, scVI-latent vs native; overall, bio and batch; per-algorithm slopes (emtrends);
#   LIGER UINMF vs LIGER (diff ~ removed + (1 + removed || tissue)); crossover of UINMF vs scVI-latent.
#   ANOVA of the three inputs: improvement ~ input * removed * algorithm + tissue random terms; Tukey pairs per % removed.
#   Run per 10% removed (removal_units).
# Section 3.2.4, number of shared genes (raw count in the dataset): one model for the three inputs on improvement over Unintegrated,
#   impr ~ input * n + input * tissue + (1 | algorithm) + (1 | algorithm:input) + (1 | dataset) + (1 | dataset:input);
#   input x n test, slope per input, slopes within each tissue as a check; overall, bio and batch. Output input_genes_model<suffix>.
# Output: figures/orthologue_input/orthologue_input<suffix>.{pdf,png,csv}, uinmf_low_shared<suffix>.{pdf,png,csv}, input_anova<suffix>.{pdf,png,csv}, uinmf_report<suffix>.{pdf,png} (report figure), input_genes_model<suffix>.{pdf,png,csv}; printout.
# Needs only CRAN stats packages (no Seurat): run with R 4.5 as
#   Rscript --vanilla analysis/orthologue_input.R
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr); library(purrr); library(stringr); library(ggplot2); library(patchwork)})

# ---------------- config ----------------
out           = 'figures/orthologue_input'
definitions   = setNames(list(list(bio = 'MAP', batch = 'SAS'), list(bio = character(0), batch = character(0))), c('_MAP_SAS', ''))   # suffix -> extra metrics
algorithms    = c('Harmony', 'FastMNN', 'RPCA', 'CCA', 'Scanorama')   # matched algorithms, as section 2.4
inputs        = c('scVI-latent', 'UINMF')
cases         = c(dissected = 0, genes2 = 20, genes4 = 40, genes6 = 60, genes8 = 80)   # % orthologues removed
reference     = 'Unintegrated'
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI')
extra_file    = 'extra_metrics/extra_metrics_all.csv'
sas_file      = 'sas/sas_all.csv'
tissue_labels = c(bone_marrow = 'Bone marrow', colon = 'Colon', breast_cancer = 'Breast cancer', fetal_liver = 'Fetal liver')
tissue_cols   = c(bone_marrow = '#E69F00', colon = '#009E73', breast_cancer = '#CC79A7', fetal_liver = '#56B4E9')   # as analysis/benchmark_figures.R
tissue_shapes = c(bone_marrow = 16, colon = 17, breast_cancer = 15, fetal_liver = 18)
# number of shared one-to-one orthologues per dataset (n_vars of the dataset's 03_concatenated_data .h5ad; 0% = dissected)
shared_genes  = tibble(tissue = rep(names(tissue_labels), each = length(cases)), case = rep(names(cases), length(tissue_labels)),
                       genes = c(123, 98, 74, 49, 25,  46, 37, 28, 18, 9,  122, 98, 73, 49, 24,  123, 98, 74, 49, 25))
removal_units = c(percent = '')   # 3.2.4 paired models / ANOVA per 10% removed (file tag); the number of shared genes has its own model below
genes_formula = impr ~ input * n + input * tissue + (1 | algorithm) + (1 | algorithm:input) + (1 | dataset) + (1 | dataset:input)   # 3.2.4 genes model
score_labels  = c(overall = 'Overall', bio = 'Bio-conservation', batch = 'Batch correction')
genes_ref     = 123                                  # 3.2.4 genes model: fitted values reported at this many shared genes (and at 25)
genes_algorithms = c('Harmony', 'BBKNN', 'FastMNN', 'RPCA', 'CCA', 'Scanorama')   # 3.2.4 genes model: as section 2.4 (BBKNN graph-only:
                                                    # scored on its own metrics, NMI / ARI / cLISI and iLISI / graph_conn, on every input)
low_genes     = 25                                   # 3.2.4 genes model: inputs compared at this many shared genes (lowest level within every tissue's range)
# -----------------------------------------
dir.create(out, showWarnings = FALSE, recursive = TRUE)
score_dir = function(tissue, case) {   # as analysis/benchmark_figures.R
  if (tissue == 'breast_cancer')
    return(sprintf('breast_cancer/data/%s/', if (case %in% c('full', 'dissected')) paste0(case, '_mouse') else case))
  if (case %in% c('full', 'dissected')) sprintf('%s/data/05_evaluated_score/%s_mouse/', tissue, case)
  else sprintf('%s/data/%s/', tissue, case)
}
read_metrics = function(f) { df = suppressMessages(read_csv(f, show_col_types = FALSE)); tibble(metric = df[[1]], value = as.numeric(df[[2]]), method = names(df)[2]) }
score_with = function(m, bio_m, batch_m) m |>   # as analysis/benchmark_figures.R
  mutate(group = case_when(metric %in% c(batch_metrics, batch_m) ~ 'batch', metric %in% c(bio_metrics, bio_m) ~ 'bio')) |>
  filter(!is.na(group)) |>
  group_by(tissue, case, method, group) |> summarise(score = mean(value), n_metrics = n(), .groups = 'drop') |>
  pivot_wider(names_from = group, values_from = c(score, n_metrics)) |>
  rename(bio = score_bio, batch = score_batch) |>
  mutate(overall = 0.6 * bio + 0.4 * batch)

met = expand_grid(tissue = names(tissue_labels), case = names(cases)) |>
  pmap_dfr(function(tissue, case) list.files(score_dir(tissue, case), pattern = '^05\\..*\\.csv$', full.names = TRUE) |>
             map_dfr(read_metrics) |> mutate(tissue = tissue, case = case)) |>
  filter(!str_detect(method, 'scIB'), !is.na(value))
extra = read_csv(extra_file, show_col_types = FALSE) |> select(tissue, case = dataset, method, MAP, NC) |>
  full_join(read_csv(sas_file, show_col_types = FALSE) |> select(tissue, case = dataset, method, SAS = sas), by = c('tissue', 'case', 'method')) |>
  mutate(case = str_remove(case, '_mouse$')) |> filter(case %in% names(cases)) |>
  pivot_longer(c(MAP, NC, SAS), names_to = 'metric', values_to = 'value') |> filter(!is.na(value))
met = bind_rows(met, extra)
methods_used = c(algorithms, as.vector(outer(algorithms, inputs, paste)))

for (suffix in names(definitions)) {
  def = definitions[[suffix]]
  cat(sprintf('\n\n################ Score definition: %s ################\n', if (suffix == '') 'reference (scib only)' else '+MAP (bio) +SAS (batch)'))
  need = c(def$bio, def$batch)
  if (length(need)) {
    have = met |> filter(metric %in% need, method %in% methods_used) |> dplyr::count(case, metric)
    full = length(methods_used) * length(tissue_labels)
    incomplete = have |> filter(n < full) |> bind_rows(anti_join(expand_grid(case = names(cases), metric = need), have, by = c('case', 'metric')))
    if (nrow(incomplete)) { cat('Skipped: extra metrics incomplete for', paste(unique(incomplete$case), collapse = ', '), '\n'); next }
  }
  sc = score_with(met, def$bio, def$batch) |> filter(method %in% methods_used)
  stopifnot(nrow(sc) == length(methods_used) * length(tissue_labels) * length(cases))
  nat = sc |> filter(method %in% algorithms) |> select(tissue, case, algorithm = method, native = overall)
  g = sc |> filter(!method %in% algorithms) |> mutate(algorithm = word(method, 1), input = str_remove(method, '^\\S+ ')) |>
    select(tissue, case, algorithm, input, overall) |> left_join(nat, by = c('tissue', 'case', 'algorithm')) |>
    mutate(gain = overall - native, removed = cases[case], r10 = removed / 10,
           algorithm = factor(algorithm, levels = algorithms), tissue = factor(tissue, levels = names(tissue_labels)), dataset = paste(tissue, case))
  write_csv(g, file.path(out, paste0('orthologue_input', suffix, '.csv')))

  fits = list(); pred = list()
  for (lv in inputs) {
    x = filter(g, input == lv)
    f = suppressMessages(lmerTest::lmer(gain ~ algorithm * r10 + (1 + r10 || tissue) + (1 | dataset), data = x,
                                        contrasts = list(algorithm = 'contr.sum')))
    cf = summary(f, ddf = 'Kenward-Roger')$coefficients
    ci = function(t) cf[t, 'Estimate'] + c(-1, 1) * qt(0.975, cf[t, 'df']) * cf[t, 'Std. Error']
    a = anova(f, type = 3, ddf = 'Kenward-Roger')
    em = suppressMessages(as.data.frame(emmeans::emmeans(f, ~ r10, at = list(r10 = c(0, 8)), lmer.df = 'kenward-roger')))
    vc = as.data.frame(lme4::VarCorr(f))
    cat(sprintf('\n== %s minus native (%d algorithms x 4 tissues x 5 cases = %d) ==\n', lv, length(algorithms), nrow(x)))
    cat(sprintf('  gain at 0%% removed: %+.4f [%+.4f, %+.4f]\n', cf['(Intercept)', 'Estimate'], ci('(Intercept)')[1], ci('(Intercept)')[2]))
    cat(sprintf('  slope per 10%% removed: %+.4f [%+.4f, %+.4f], df %.1f, p = %.3g\n', cf['r10', 'Estimate'], ci('r10')[1], ci('r10')[2],
                cf['r10', 'df'], cf['r10', 'Pr(>|t|)']))
    cat(sprintf('  predicted gain at 0%% / 80%%: %+.4f [%+.4f, %+.4f] / %+.4f [%+.4f, %+.4f]\n', em$emmean[1], em$lower.CL[1], em$upper.CL[1],
                em$emmean[2], em$lower.CL[2], em$upper.CL[2]))
    cat(sprintf('  algorithm x removal: F(%d, %.1f) = %.2f, p = %.3g (do algorithms lose the advantage at different rates?)\n',
                a['algorithm:r10', 'NumDF'], a['algorithm:r10', 'DenDF'], a['algorithm:r10', 'F value'], a['algorithm:r10', 'Pr(>F)']))
    cat('  random-effect SDs: ', paste(sprintf('%s %.4f', paste(vc$grp, ifelse(is.na(vc$var1), '', vc$var1)), vc$sdcor), collapse = ', '),
        if (lme4::isSingular(f)) '(singular)' else '', '\n')
    # curvature check: removal as a 5-level factor vs linear (ML fits)
    f_lin = suppressMessages(lme4::lmer(gain ~ algorithm * r10 + (1 + r10 || tissue) + (1 | dataset), data = x, REML = FALSE))
    f_fac = suppressMessages(lme4::lmer(gain ~ algorithm * r10 + factor(removed) + (1 + r10 || tissue) + (1 | dataset), data = x, REML = FALSE))
    lr = anova(f_lin, f_fac)
    cat(sprintf('  departure from linear trend: chi2(%d) = %.2f, p = %.3f\n', lr$Df[2], lr$Chisq[2], lr$`Pr(>Chisq)`[2]))
    per = x |> group_by(tissue) |> summarise(slope = coef(lm(gain ~ algorithm + r10))['r10'], at0 = mean(gain[removed == 0]), at80 = mean(gain[removed == 80]))
    cat('  per tissue (slope per 10%, mean gain at 0% and 80%):\n'); print(as.data.frame(per |> mutate(across(where(is.double), \(v) round(v, 4)))), row.names = FALSE)
    pa = x |> group_by(algorithm) |> summarise(slope = coef(lm(gain ~ tissue + r10))['r10'], at0 = mean(gain[removed == 0]), at80 = mean(gain[removed == 80]))
    cat('  per algorithm:\n'); print(as.data.frame(pa |> mutate(across(where(is.double), \(v) round(v, 4)))), row.names = FALSE)
    cat(sprintf('  gain > 0 in %d / %d algorithm x tissue x case\n', sum(x$gain > 0), nrow(x)))
    pred[[lv]] = suppressMessages(as.data.frame(emmeans::emmeans(f, ~ r10, at = list(r10 = seq(0, 8, 0.5)), lmer.df = 'kenward-roger'))) |> mutate(input = lv)
  }

  # ---- figure: mean gain over the algorithms per tissue and case; fitted mean with 95% CI ----
  tm = g |> group_by(input, tissue, removed) |> summarise(gain = mean(gain), .groups = 'drop')
  pr = bind_rows(pred) |> mutate(removed = r10 * 10)
  lab_in = c(`scVI-latent` = 'scVI-latent minus native', UINMF = 'UINMF minus native')
  p = ggplot(tm, aes(removed, gain)) +
    geom_hline(yintercept = 0, colour = 'grey40', linetype = 'dashed') +
    geom_ribbon(data = pr, aes(y = emmean, ymin = lower.CL, ymax = upper.CL, fill = 'fit'), alpha = 0.7) +
    geom_line(data = pr, aes(y = emmean, linetype = 'fit'), linewidth = 1.1) +
    geom_line(aes(colour = tissue), linewidth = 0.6, alpha = 0.8) +
    geom_point(aes(colour = tissue, shape = tissue), size = 2.4) +
    facet_wrap(~ input, labeller = as_labeller(lab_in)) +
    scale_colour_manual(values = tissue_cols, labels = tissue_labels, name = NULL) +
    scale_shape_manual(values = tissue_shapes, labels = tissue_labels, name = NULL) +
    scale_fill_manual(values = c(fit = 'grey80'), labels = c(fit = 'Fitted mean over algorithms\n(95% CI)'), name = NULL) +
    scale_linetype_manual(values = c(fit = 'solid'), labels = c(fit = 'Fitted mean over algorithms\n(95% CI)'), name = NULL) +
    scale_x_continuous(breaks = cases, labels = \(v) paste0(v, '%')) +
    labs(x = 'Shared orthologues removed (0% = dissected)', y = 'Gain over the same algorithm\non native input (overall score)') +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = 'grey95', colour = 'grey70'),
          strip.text = element_text(face = 'bold'), legend.position = 'bottom')
  for (ext in c('pdf', 'png'))
    ggsave(file.path(out, paste0('orthologue_input', suffix, '.', ext)), p, width = 10, height = 5, dpi = 300, bg = 'white',
           device = if (ext == 'pdf') cairo_pdf else 'png')

  # ===================== 3.2.4 Is UINMF the answer to low shared genes? =====================
  # Removal is measured in each unit of `removal_units`: 'percent' = per 10% of shared orthologues removed; 'genes' = the
  # number of shared genes in the dataset (raw count, slope per shared gene), the same scale for every tissue, so colon's
  # smaller panel (46 to 9 genes) enters through its gene counts. Fitted values in 'genes' at `genes_ref` and 25 shared genes.
  ref = score_with(met, def$bio, def$batch) |> filter(method == reference) |> select(tissue, case, ref = overall)
  stopifnot(nrow(ref) == length(tissue_labels) * length(cases))
  refs = score_with(met, def$bio, def$batch) |> filter(method == reference) |> select(tissue, case, overall, bio, batch) |>
    pivot_longer(c(overall, bio, batch), names_to = 'score', values_to = 'ref')
  g0 = shared_genes |> filter(case == 'dissected') |> select(tissue, genes0 = genes)
  r3 = \(v) round(v, 4)
  in_cols = c(Native = '#4477AA', `scVI-latent` = '#7B5EA7', UINMF = '#CCBB44')
  for (unit in names(removal_units)) {
    tag = removal_units[[unit]]
    genes_unit = unit == 'genes'
    cat(sprintf('\n\n======== 3.2.4 Is UINMF the answer to low shared genes? (removal %s) ========\n',
                if (genes_unit) 'as the number of shared genes in the dataset' else 'per 10% of shared orthologues'))
    xv = shared_genes |> left_join(g0, by = 'tissue') |>
      # model predictor ('r10' kept as the name in both units): percent = per 10% removed (0 at the dissected data);
      # genes = number of shared genes in the dataset (raw count)
      mutate(r10 = if (genes_unit) genes else cases[case] / 10, base = 0) |> select(tissue, case, genes, r10, base)
    base_term = ''
    base_term_av = ''
    at_ref = list()
    r_0 = if (genes_unit) genes_ref else 0                                      # reported start: 123 shared genes / 0% removed
    r_end = if (genes_unit) round(0.2 * genes_ref) else 8                       # reported end: 25 shared genes / 80% removed
    r_curve = if (genes_unit) min(shared_genes$genes) else r_end                # fitted lines drawn to here (genes: colon's 9)
    xpos = \(r) if (genes_unit) r else 10 * r                                   # model scale -> plot axis
    per = if (genes_unit) 'per shared gene' else 'per 10% removed'
    fmt = if (genes_unit) '%+.5f' else '%+.4f'; r3 = \(v) round(v, if (genes_unit) 5 else 4)
    at_lab = if (genes_unit) sprintf('%d and %d shared genes', genes_ref, round(xpos(r_end))) else '0% and 80% removed'
    x_scale = if (genes_unit) scale_x_reverse(breaks = seq(0, 125, 25)) else scale_x_continuous(breaks = cases, labels = \(v) paste0(v, '%'))
    x_lab = if (genes_unit) 'Number of shared orthologues' else 'Shared orthologues removed'

    wide = sc |> mutate(algorithm = word(method, 1), input = ifelse(method %in% algorithms, 'Native', str_remove(method, '^\\S+ '))) |>
      select(tissue, case, algorithm, input, overall, bio, batch) |>
      pivot_longer(c(overall, bio, batch), names_to = 'score') |> pivot_wider(names_from = input, values_from = value) |>
      left_join(xv, by = c('tissue', 'case')) |>
      mutate(removed = cases[case], algorithm = factor(algorithm, levels = algorithms),
             tissue = factor(tissue, levels = names(tissue_labels)), dataset = paste(tissue, case))
    stopifnot(!anyNA(wide$r10))
    comps = list(`UINMF vs native` = c('UINMF', 'Native'), `UINMF vs scVI-latent` = c('UINMF', 'scVI-latent'),
                 `scVI-latent vs native` = c('scVI-latent', 'Native'))
    fit_c = function(x) suppressMessages(lmerTest::lmer(as.formula(paste0('diff ~ algorithm * r10', base_term, ' + (1 + r10 || tissue) + (1 | dataset)')),
                                                        data = x, contrasts = list(algorithm = 'contr.sum')))
    summ = function(f) {
      cf = summary(f, ddf = 'Kenward-Roger')$coefficients; a = anova(f, type = 3, ddf = 'Kenward-Roger')
      em = suppressMessages(as.data.frame(emmeans::emmeans(f, ~ r10, at = c(list(r10 = c(r_0, r_end)), at_ref), lmer.df = 'kenward-roger')))
      sl = cf['r10', 'Estimate'] + c(-1, 1) * qt(0.975, cf['r10', 'df']) * cf['r10', 'Std. Error']
      tibble(at0 = em$emmean[1], at0_lo = em$lower.CL[1], at0_hi = em$upper.CL[1], slope = cf['r10', 'Estimate'], slope_lo = sl[1], slope_hi = sl[2],
             p_slope = cf['r10', 'Pr(>|t|)'], at80 = em$emmean[2], at80_lo = em$lower.CL[2], at80_hi = em$upper.CL[2],
             p_alg_slope = a['algorithm:r10', 'Pr(>F)'], p_base = NA_real_)
    }
    res = list(); trends = list(); preds = list()
    for (cn in names(comps)) for (s_ in c('overall', 'bio', 'batch')) {
      x = wide |> filter(score == s_) |> mutate(diff = .data[[comps[[cn]][1]]] - .data[[comps[[cn]][2]]])
      f = fit_c(x)
      res[[paste(cn, s_)]] = summ(f) |> mutate(comparison = cn, score = s_, .before = 1)
      trends[[paste(cn, s_)]] = suppressMessages(as.data.frame(emmeans::emtrends(f, ~ algorithm, var = 'r10', lmer.df = 'kenward-roger'))) |>
        mutate(comparison = cn, score = s_)
      if (s_ == 'overall') preds[[cn]] = suppressMessages(as.data.frame(emmeans::emmeans(f, ~ r10, at = c(list(r10 = seq(r_0, r_curve, length.out = 17)), at_ref),
                                                                                          lmer.df = 'kenward-roger'))) |> mutate(comparison = cn)
    }
    res = bind_rows(res); trends = bind_rows(trends)
    # LIGER: the algorithm UINMF was designed for (native LIGER vs LIGER UINMF), its own model (one algorithm)
    lig = score_with(met, def$bio, def$batch) |> filter(method %in% c('LIGER', 'LIGER UINMF')) |>
      select(tissue, case, method, overall, bio, batch) |> pivot_longer(c(overall, bio, batch), names_to = 'score') |>
      pivot_wider(names_from = method, values_from = value) |> left_join(xv, by = c('tissue', 'case')) |>
      mutate(diff = `LIGER UINMF` - LIGER, tissue = factor(tissue, levels = names(tissue_labels)))
    stopifnot(nrow(lig) == 3 * 4 * length(cases))
    lig_res = bind_rows(lapply(c('overall', 'bio', 'batch'), \(s_) {
      f = suppressMessages(lmerTest::lmer(as.formula(paste0('diff ~ r10', base_term, ' + (1 + r10 || tissue)')), data = filter(lig, score == s_)))
      cf = summary(f, ddf = 'Kenward-Roger')$coefficients
      em = suppressMessages(as.data.frame(emmeans::emmeans(f, ~ r10, at = c(list(r10 = r_0), at_ref), lmer.df = 'kenward-roger')))
      tibble(score = s_, at0 = em$emmean, slope = cf['r10', 1], slope_lo = cf['r10', 1] - qt(0.975, cf['r10', 'df']) * cf['r10', 2],
             slope_hi = cf['r10', 1] + qt(0.975, cf['r10', 'df']) * cf['r10', 2], p_slope = cf['r10', 'Pr(>|t|)'],
             pos = sum(filter(lig, score == s_)$diff > 0))
    }))
    write_csv(bind_rows(res, lig_res |> mutate(comparison = 'LIGER UINMF vs LIGER')), file.path(out, paste0('uinmf_low_shared', tag, suffix, '.csv')))
    cat(sprintf('Mean over the matched algorithms (fitted difference [95%% CI] at %s; slope %s):\n', at_lab, per))
    res |> transmute(comparison, score, at_start = sprintf('%+.3f [%+.3f, %+.3f]', at0, at0_lo, at0_hi),
                     slope = sprintf(paste0(fmt, ' [', fmt, ', ', fmt, ']'), slope, slope_lo, slope_hi), p_slope = signif(p_slope, 2),
                     at_end = sprintf('%+.3f [%+.3f, %+.3f]', at80, at80_lo, at80_hi), p_alg_x_removal = signif(p_alg_slope, 2),
                     p_tissue_start_genes = signif(p_base, 2)) |>
      select(where(\(v) !all(is.na(v)))) |> as.data.frame() |> print(row.names = FALSE)
    uv = filter(res, comparison == 'UINMF vs scVI-latent', score == 'overall')
    cross = r_0 - uv$at0 / uv$slope   # model scale
    cross_ok = cross >= min(r_0, r_curve) && cross <= max(r_0, r_curve)
    cross_txt = if (genes_unit) sprintf('%.0f shared genes', xpos(cross)) else sprintf('%.0f%% removed', xpos(cross))
    cat(sprintf('UINMF vs scVI-latent (overall): fitted difference crosses 0 at %s\n', if (cross_ok) cross_txt else 'no point in the range'))
    cat(sprintf('\nSlope %s, by algorithm (overall; emtrends, 95%% CI):\n', per))
    trends |> filter(score == 'overall', comparison != 'scVI-latent vs native') |>
      transmute(comparison, algorithm, slope = r3(r10.trend), lo = r3(lower.CL), hi = r3(upper.CL)) |> as.data.frame() |> print(row.names = FALSE)
    cat(sprintf('\nLIGER UINMF vs LIGER (diff ~ removed%s + (1 + removed || tissue)):\n', base_term))
    lig_res |> mutate(across(where(is.double), \(v) signif(v, 3))) |> as.data.frame() |> print(row.names = FALSE)
    if (!genes_unit) {
      cat('\nRaw mean differences over algorithms and tissues, by % removed (overall):\n')
      wide |> filter(score == 'overall') |> group_by(removed) |>
        summarise(`UINMF - native` = r3(mean(UINMF - Native)), `UINMF - scVI-latent` = r3(mean(UINMF - `scVI-latent`)),
                  `scVI-latent - native` = r3(mean(`scVI-latent` - Native)),
                  `UINMF best of 3 inputs` = sprintf('%d/20', sum(UINMF > pmax(Native, `scVI-latent`)))) |> as.data.frame() |> print(row.names = FALSE)
    }

    # ---- 3.2.4 ANOVA: do the three inputs differ, and does that change as orthologues are removed? ----
    # outcome = improvement over Unintegrated in the same tissue and case; input, removal (5 levels) and algorithm fixed
    # (full factorial, sum-to-zero), tissue random with its interactions: input x tissue, algorithm x tissue, dataset
    # (= tissue:case) and input x dataset (error term for input x removal). Type III, Kenward-Roger.
    # The 5-level factor model does not depend on the unit: printed and plotted once (percent); the linear version per unit.
    av = sc |> mutate(algorithm = word(method, 1), input = ifelse(method %in% algorithms, 'Native', str_remove(method, '^\\S+ '))) |>
      left_join(ref, by = c('tissue', 'case')) |> left_join(xv, by = c('tissue', 'case')) |>
      mutate(impr = overall - ref, input = factor(input, levels = c('Native', inputs)), algorithm = factor(algorithm, levels = algorithms),
             tissue = factor(tissue, levels = names(tissue_labels)), removed = factor(cases[case]), dataset = paste(tissue, case))
    stopifnot(nrow(av) == length(algorithms) * 3 * length(tissue_labels) * length(cases), !anyNA(av$impr), !anyNA(av$r10))
    rand = '(1 | tissue) + (1 | tissue:input) + (1 | tissue:algorithm) + (1 | dataset) + (1 | dataset:input)'
    ctr = list(input = 'contr.sum', algorithm = 'contr.sum', removed = 'contr.sum')
    if (!genes_unit) {
      fa = suppressMessages(lmerTest::lmer(as.formula(paste('impr ~ input * removed * algorithm +', rand)), data = av, contrasts = ctr))
      at = anova(fa, type = 3, ddf = 'Kenward-Roger')
      cat(sprintf('\nANOVA: improvement over Unintegrated ~ input * removed * algorithm + %s\n', rand))
      cat(sprintf('(%d algorithms x 3 inputs x 4 tissues x 5 cases = %d; removal as a 5-level factor; type III, Kenward-Roger)\n', length(algorithms), nrow(av)))
      print(as.data.frame(at) |> transmute(NumDF, DenDF = round(DenDF, 1), F = round(`F value`, 2), p = signif(`Pr(>F)`, 3)))
      vc = as.data.frame(lme4::VarCorr(fa))
      cat('random-effect SDs: ', paste(sprintf('%s %.4f', vc$grp, vc$sdcor), collapse = ', '), if (lme4::isSingular(fa)) '(singular)' else '', '\n')
    }
    # linear version of input x removal (2 df instead of 8)
    fl = suppressMessages(lmerTest::lmer(as.formula(paste0('impr ~ input * r10 * algorithm', base_term_av, ' + ', rand)), data = av,
                                         contrasts = ctr[c('input', 'algorithm')]))
    al = anova(fl, type = 3, ddf = 'Kenward-Roger')
    cat(sprintf('\nremoval as linear (%s%s): input x removal F(%d, %.1f) = %.2f, p = %.3g%s\n', per,
                '',
                al['input:r10', 'NumDF'], al['input:r10', 'DenDF'], al['input:r10', 'F value'], al['input:r10', 'Pr(>F)'],
                ''))
    if (!genes_unit) {
      emi = suppressMessages(emmeans::emmeans(fa, ~ input, lmer.df = 'kenward-roger'))
      cat('Input means over all removal levels and algorithms, with Tukey-adjusted pairwise differences:\n')
      print(as.data.frame(emi) |> transmute(input, mean = r3(emmean), lo = r3(lower.CL), hi = r3(upper.CL)), row.names = FALSE)
      print(as.data.frame(pairs(emi, adjust = 'tukey')) |> transmute(contrast, diff = r3(estimate), SE = r3(SE), df = round(df, 1), p = signif(p.value, 3)), row.names = FALSE)
      emr = suppressMessages(emmeans::emmeans(fa, ~ input | removed, lmer.df = 'kenward-roger'))
      pr_r = as.data.frame(pairs(emr, adjust = 'tukey')) |> mutate(contrast = gsub('[()]', '', contrast))
      ci_r = as.data.frame(confint(pairs(emr, adjust = 'tukey'))) |> mutate(contrast = gsub('[()]', '', contrast))
      cat('Pairwise differences between inputs at each % removed (Tukey within each level):\n')
      pr_r |> left_join(ci_r |> select(contrast, removed, lower.CL, upper.CL), by = c('contrast', 'removed')) |>
        transmute(removed, contrast, diff = sprintf('%+.3f [%+.3f, %+.3f]', estimate, lower.CL, upper.CL), df = round(df, 1), p = signif(p.value, 3)) |>
        print(row.names = FALSE)
      write_csv(bind_rows(as.data.frame(at) |> tibble::rownames_to_column('term') |> mutate(table = 'anova'),
                          pr_r |> left_join(ci_r |> select(contrast, removed, lower.CL, upper.CL), by = c('contrast', 'removed')) |>
                            mutate(table = 'pairwise_by_removed', removed = as.character(removed))),
                file.path(out, paste0('input_anova', suffix, '.csv')))

      # ---- ANOVA figure: (A) input means by % removed; (B) pairwise differences by % removed ----
      ea = as.data.frame(emr) |> mutate(removed = as.numeric(as.character(removed)))
      raw = av |> group_by(input, tissue, removed) |> summarise(impr = mean(impr), .groups = 'drop') |> mutate(removed = as.numeric(as.character(removed)))
      pa1 = ggplot(ea, aes(removed, emmean, colour = input)) +
        geom_hline(yintercept = 0, colour = 'grey40', linetype = 'dashed') +
        geom_point(data = raw, aes(y = impr, shape = tissue), position = position_dodge(width = 9), size = 1.4, alpha = 0.45) +
        geom_line(position = position_dodge(width = 9), linewidth = 0.7) +
        geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 9), size = 0.5, linewidth = 0.8) +
        scale_colour_manual(values = in_cols, name = 'Input') +
        scale_shape_manual(values = tissue_shapes, labels = tissue_labels, name = 'Tissue (mean over algorithms)') +
        scale_x_continuous(breaks = cases, labels = \(v) paste0(v, '%')) +
        labs(x = 'Shared orthologues removed', y = 'Improvement over Unintegrated\n(estimated mean over algorithms, 95% CI)') +
        theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = 'bottom', legend.box = 'vertical')
      pb = pr_r |> left_join(ci_r |> select(contrast, removed, lower.CL, upper.CL), by = c('contrast', 'removed')) |>
        mutate(removed = factor(paste0(removed, '%'), levels = rev(paste0(cases, '%'))), contrast = factor(contrast, levels = unique(contrast)))
      pa2 = ggplot(pb, aes(estimate, removed)) +
        geom_vline(xintercept = 0, colour = 'grey40') +
        geom_pointrange(aes(xmin = lower.CL, xmax = upper.CL, colour = p.value < 0.05), size = 0.4) +
        facet_wrap(~ contrast, ncol = 1) +
        scale_colour_manual(values = c(`TRUE` = 'black', `FALSE` = 'grey60'), labels = c(`TRUE` = 'p < 0.05', `FALSE` = 'p ≥ 0.05'), name = 'Tukey-adjusted') +
        labs(x = 'Difference in improvement\n(Tukey-adjusted 95% CI)', y = 'Orthologues removed') +
        theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = 'grey95', colour = 'grey70'),
                                         strip.text = element_text(face = 'bold'), legend.position = 'bottom')
      pAN = (pa1 | pa2) + patchwork::plot_layout(widths = c(1.5, 1)) + patchwork::plot_annotation(tag_levels = 'A') &
        theme(plot.tag = element_text(face = 'bold', size = 13))
      for (ext in c('pdf', 'png'))
        ggsave(file.path(out, paste0('input_anova', suffix, '.', ext)), pAN, width = 12, height = 6.5, dpi = 300, bg = 'white',
               device = if (ext == 'pdf') cairo_pdf else 'png')
    }

    # ---- figure: (A) UINMF minus scVI-latent by tissue + fit; (B) slope by algorithm; (C) slope by score ----
    tm2 = wide |> filter(score == 'overall') |> group_by(tissue, removed, genes) |> summarise(diff = mean(UINMF - `scVI-latent`), .groups = 'drop') |>
      mutate(x = if (genes_unit) genes else removed)
    pr2 = preds[['UINMF vs scVI-latent']] |> mutate(x = xpos(r10))
    fit_lab = 'Fitted mean over algorithms (95% CI)'
    pA = ggplot(tm2, aes(x, diff)) +
      geom_hline(yintercept = 0, colour = 'grey40', linetype = 'dashed') +
      geom_ribbon(data = pr2, aes(y = emmean, ymin = lower.CL, ymax = upper.CL, fill = 'fit'), alpha = 0.7) +
      geom_line(data = pr2, aes(y = emmean, linetype = 'fit'), linewidth = 1.1) +
      geom_line(aes(colour = tissue), linewidth = 0.6, alpha = 0.8) + geom_point(aes(colour = tissue, shape = tissue), size = 2.4) +
      scale_colour_manual(values = tissue_cols, labels = tissue_labels, name = NULL) +
      scale_shape_manual(values = tissue_shapes, labels = tissue_labels, name = NULL) +
      scale_fill_manual(values = c(fit = 'grey80'), labels = c(fit = fit_lab), name = NULL) +
      scale_linetype_manual(values = c(fit = 'solid'), labels = c(fit = fit_lab), name = NULL) +
      x_scale + labs(x = x_lab, y = 'UINMF minus scVI-latent input\n(overall score, same algorithm)') +
      theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = 'bottom', legend.box = 'vertical')
    comp_cols = c(`UINMF vs native` = unname(in_cols['Native']), `UINMF vs scVI-latent` = unname(in_cols['scVI-latent']))   # colour = the input UINMF is compared with (as the report figure)
    tb = trends |> filter(score == 'overall', comparison %in% names(comp_cols)) |>
      transmute(comparison, row = as.character(algorithm), slope = r10.trend, lo = lower.CL, hi = upper.CL) |>
      bind_rows(lig_res |> filter(score == 'overall') |> transmute(comparison = 'UINMF vs native', row = 'LIGER', slope, lo = slope_lo, hi = slope_hi)) |>
      mutate(row = factor(row, levels = rev(c(algorithms, 'LIGER')), labels = rev(c(ifelse(algorithms %in% c('CCA', 'RPCA'), paste0(algorithms, ' (Seurat)'), algorithms), 'LIGER'))))
    pB = ggplot(tb, aes(slope, row, colour = comparison)) +
      geom_vline(xintercept = 0, colour = 'grey40') +
      geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.55), size = 0.45) +
      scale_colour_manual(values = comp_cols, name = NULL) +
      labs(x = sprintf('Change in difference %s\n(overall score, 95%% CI)', per), y = NULL) +
      theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = 'bottom', legend.direction = 'vertical',
                                       plot.margin = margin(5.5, 18, 5.5, 5.5))
    tc = res |> filter(comparison %in% names(comp_cols)) |>
      mutate(score = factor(score, levels = rev(c('overall', 'bio', 'batch')), labels = rev(c('Overall', 'Bio-conservation', 'Batch correction'))))
    pC = ggplot(tc, aes(slope, score, colour = comparison)) +
      geom_vline(xintercept = 0, colour = 'grey40') +
      geom_pointrange(aes(xmin = slope_lo, xmax = slope_hi), position = position_dodge(width = 0.55), size = 0.45) +
      scale_colour_manual(values = comp_cols, name = NULL, guide = 'none') +
      labs(x = sprintf('Change %s\n(mean over algorithms, 95%% CI)', per), y = NULL) +
      theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
    pU = (pA | (pB / pC + patchwork::plot_layout(heights = c(2, 1)))) + patchwork::plot_layout(widths = c(1.25, 1)) +
      patchwork::plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 13))
    for (ext in c('pdf', 'png'))
      ggsave(file.path(out, paste0('uinmf_low_shared', tag, suffix, '.', ext)), pU, width = 12.5, height = 6.5, dpi = 300, bg = 'white',
             device = if (ext == 'pdf') cairo_pdf else 'png')

    # ---- 3.2.4 report figure: (A) input means (percent: 5-level ANOVA; genes: linear model) ; (B) slope per input for
    # overall / bio / batch (linear model, each score's improvement over Unintegrated); (C) UINMF minus scVI-latent +
    # crossover; (D) by algorithm ----
    avs = av |> select(tissue, case, algorithm, input, removed, r10, base, dataset, method) |>
      left_join(sc |> select(tissue, case, method, overall, bio, batch), by = c('tissue', 'case', 'method')) |>
      pivot_longer(c(overall, bio, batch), names_to = 'score') |> left_join(refs, by = c('tissue', 'case', 'score')) |> mutate(impr = value - ref)
    slopes = map_dfr(c('overall', 'bio', 'batch'), \(s_) {
      f = suppressMessages(lmerTest::lmer(as.formula(paste0('impr ~ input * r10 * algorithm', base_term_av, ' + ', rand)), data = filter(avs, score == s_),
                                          contrasts = ctr[c('input', 'algorithm')]))
      a = anova(f, type = 3, ddf = 'Kenward-Roger')
      suppressMessages(as.data.frame(emmeans::emtrends(f, ~ input, var = 'r10', lmer.df = 'kenward-roger'))) |>
        mutate(score = s_, p_int = a['input:r10', 'Pr(>F)'])
    })
    cat(sprintf('\nSlope %s by input and score (improvement over Unintegrated; input x removal p from the linear model):\n', per))
    slopes |> transmute(score, input, slope = r3(r10.trend), lo = r3(lower.CL), hi = r3(upper.CL), p_input_x_removal = signif(p_int, 3)) |>
      print(row.names = FALSE)
    if (genes_unit) {
      ea = suppressMessages(as.data.frame(emmeans::emmeans(fl, ~ input | r10, at = c(list(r10 = seq(r_0, r_curve, length.out = 17)), at_ref),
                                                           lmer.df = 'kenward-roger'))) |> mutate(x = xpos(r10))
      ann = sprintf('Input × shared genes: F(%d, %.0f) = %.2f, p = %.3f', al['input:r10', 'NumDF'], al['input:r10', 'DenDF'],
                    al['input:r10', 'F value'], al['input:r10', 'Pr(>F)'])
      rA = ggplot(ea, aes(x, emmean, colour = input, fill = input)) +
        geom_hline(yintercept = 0, colour = 'grey40', linetype = 'dashed') +
        geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.12, colour = NA) +
        geom_line(linewidth = 1) +
        annotate('text', x = genes_ref, y = 0.01, label = ann, hjust = 0, vjust = 0, size = 3.3) +
        scale_fill_manual(values = in_cols, name = 'Input') +
        labs(y = 'Improvement over Unintegrated\n(fitted mean over algorithms, 95% CI)')
    } else {
      ea = as.data.frame(emr) |> mutate(x = as.numeric(as.character(removed)))
      ann = sprintf('Input: F(%d, %.0f) = %.2f, p = %.2f\nInput × removal: F(%d, %.0f) = %.2f, p = %.3f',
                    at['input', 'NumDF'], at['input', 'DenDF'], at['input', 'F value'], at['input', 'Pr(>F)'],
                    at['input:removed', 'NumDF'], at['input:removed', 'DenDF'], at['input:removed', 'F value'], at['input:removed', 'Pr(>F)'])
      rA = ggplot(ea, aes(x, emmean, colour = input)) +
        geom_hline(yintercept = 0, colour = 'grey40', linetype = 'dashed') +
        geom_line(position = position_dodge(width = 5), linewidth = 0.8) +
        geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), position = position_dodge(width = 5), size = 0.5, linewidth = 0.8) +
        annotate('text', x = -4, y = 0.01, label = ann, hjust = 0, vjust = 0, size = 3.3, lineheight = 0.95) +
        labs(y = 'Improvement over Unintegrated\n(estimated mean over algorithms, 95% CI)')
    }
    rA = rA + scale_colour_manual(values = in_cols, name = 'Input') + x_scale + labs(x = x_lab) +
      theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = 'bottom')
    rB = slopes |> mutate(score = factor(score, levels = rev(c('overall', 'bio', 'batch')), labels = rev(c('Overall', 'Bio-conservation', 'Batch correction'))),
                          input = factor(input, levels = rev(levels(av$input)))) |>
      ggplot(aes(r10.trend, score, colour = input)) +
      geom_vline(xintercept = 0, colour = 'grey40') +
      geom_pointrange(aes(xmin = lower.CL, xmax = upper.CL), position = position_dodge(width = 0.6), size = 0.45) +
      scale_colour_manual(values = in_cols, guide = 'none') +
      labs(x = sprintf('Change in improvement %s\n(mean over algorithms, 95%% CI)', per), y = NULL) +
      theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
    rC = pA
    if (cross_ok) rC = rC + annotate('segment', x = xpos(cross), xend = xpos(cross), y = -Inf, yend = 0, linetype = 'dotted', colour = 'black') +
      annotate('label', x = xpos(cross), y = -0.13, label = paste0('fitted crossover\n≈ ', sub(' remaining| removed', '', cross_txt)),
               size = 3, label.size = 0, fill = 'white', lineheight = 0.9)
    comp_cols2 = c(`UINMF vs native` = unname(in_cols['Native']), `UINMF vs scVI-latent` = unname(in_cols['scVI-latent']))
    rD = suppressMessages(pB + scale_colour_manual(values = comp_cols2, name = NULL) +
      labs(x = sprintf('Change in UINMF advantage %s\n(overall score, 95%% CI)', per)) + theme(plot.margin = margin(5.5, 18, 5.5, 5.5)))
    pR = (rA + rB + patchwork::plot_layout(widths = c(1.5, 1))) / (rC + rD + patchwork::plot_layout(widths = c(1.5, 1))) +
      patchwork::plot_layout(heights = c(1, 1.15)) + patchwork::plot_annotation(tag_levels = 'A') &
      theme(plot.tag = element_text(face = 'bold', size = 13))
    for (ext in c('pdf', 'png'))
      ggsave(file.path(out, paste0('uinmf_report', tag, suffix, '.', ext)), pR, width = 13, height = 11, dpi = 300, bg = 'white',
             device = if (ext == 'pdf') cairo_pdf else 'png')
  }

  # ===================== 3.2.4 Number of shared genes: one model for the three inputs =====================
  # Improvement over Unintegrated of algorithm m on input i in dataset k (tissue j); n = shared one-to-one orthologues in the dataset:
  #   impr ~ input * n + input * tissue + (1 | algorithm) + (1 | algorithm:input) + (1 | dataset) + (1 | dataset:input)
  # Input, n and tissue fixed (sum-to-zero; tissue fixed as section 2.4), so the slope on n is estimated within tissues and each
  # tissue's starting number of shared genes is absorbed by tissue. Algorithm random (as section 2.4). dataset:input is the error
  # term for input x n (n varies only between datasets). Type III F tests, Kenward-Roger df; slope per input (emtrends) and slope
  # differences; fitted improvement averaged over tissues (emmeans). Fitted to the overall, bio and batch scores.
  # Check of the common slope: slope per input within each tissue,
  #   impr ~ input * n + algorithm * input + (1 | dataset) + (1 | dataset:input), one tissue at a time.
  # Output: input_genes_model<suffix>.{pdf,png,csv} (A fitted lines by input, B differences between inputs along n); input_genes_model_tissue_slopes<suffix> (appendix)
  methods_g = c(genes_algorithms, as.vector(outer(genes_algorithms, inputs, paste)))
  scg = score_with(met, def$bio, def$batch) |> filter(method %in% methods_g)   # metrics a method lacks (BBKNN: MAP, SAS, ...) are skipped
  stopifnot(nrow(scg) == length(methods_g) * length(tissue_labels) * length(cases))
  gd = scg |> mutate(algorithm = word(method, 1), input = ifelse(method %in% genes_algorithms, 'Native', str_remove(method, '^\\S+ '))) |>
    select(tissue, case, algorithm, input, overall, bio, batch) |> pivot_longer(c(overall, bio, batch), names_to = 'score') |>
    left_join(refs, by = c('tissue', 'case', 'score')) |> left_join(shared_genes, by = c('tissue', 'case')) |>
    mutate(impr = value - ref, n = genes, input = factor(input, levels = c('Native', inputs)), algorithm = factor(algorithm, levels = genes_algorithms),
           tissue = factor(tissue, levels = names(tissue_labels)), dataset = paste(tissue, case))
  stopifnot(nrow(gd) == 3 * 3 * length(genes_algorithms) * length(tissue_labels) * length(cases), !anyNA(gd$impr), !anyNA(gd$n))
  cat(sprintf('\n\n======== 3.2.4 Number of shared genes: one model for the three inputs (%d algorithms x 3 inputs x %d datasets) ========\n',
              length(genes_algorithms), n_distinct(gd$dataset)))
  cat(sprintf('%s\n', deparse1(genes_formula)))
  # check: on the dissected data (0% removed), the improvements equal those used in section 2.4 (analysis/benchmark_figures.R output)
  b24 = read_csv(sprintf('figures/benchmark_summary/benchmark_improvement%s.csv', suffix), show_col_types = FALSE)
  chk = gd |> filter(case == 'dissected', score == 'overall') |>
    mutate(method = ifelse(input == 'Native', as.character(algorithm), paste(algorithm, input))) |>
    inner_join(b24 |> transmute(tissue, method = str_remove(method, ' [†‡]+$'), imp_overall), by = c('tissue', 'method'))
  cat(sprintf('check vs section 2.4 input (dissected, overall): %d of %d matched, max |difference| = %.2g\n',
              nrow(chk), nrow(filter(gd, case == 'dissected', score == 'overall')), max(abs(chk$impr - chk$imp_overall))))
  g_anova = list(); g_slope = list(); g_pair = list(); g_tissue = list(); g_fit = list(); g_low = list(); g_dline = list()
  for (s_ in names(score_labels)) {
    x = filter(gd, score == s_)
    f = suppressMessages(lmerTest::lmer(genes_formula, data = x, contrasts = list(input = 'contr.sum', tissue = 'contr.sum')))
    a = anova(f, type = 3, ddf = 'Kenward-Roger')
    g_anova[[s_]] = as.data.frame(a) |> tibble::rownames_to_column('term') |> transmute(score = s_, term, NumDF, DenDF, F = `F value`, p = `Pr(>F)`)
    tr = emmeans::emtrends(f, ~ input, var = 'n', lmer.df = 'kenward-roger')
    g_slope[[s_]] = as.data.frame(summary(tr, infer = TRUE)) |> transmute(score = s_, input, slope = n.trend, lo = lower.CL, hi = upper.CL, df, p = p.value)
    g_pair[[s_]] = as.data.frame(summary(pairs(tr, adjust = 'none'), infer = TRUE)) |>
      transmute(score = s_, contrast = gsub('[()]', '', contrast), diff = estimate, lo = lower.CL, hi = upper.CL, df, p = p.value)
    g_fit[[s_]] = as.data.frame(emmeans::emmeans(f, ~ input | n, at = list(n = seq(max(x$n), min(x$n), length.out = 40)), lmer.df = 'kenward-roger')) |>
      mutate(score = s_)
    vc = as.data.frame(lme4::VarCorr(f))
    cat(sprintf('\n-- %s --\n', score_labels[[s_]]))
    print(g_anova[[s_]] |> transmute(term, NumDF, DenDF = round(DenDF, 1), F = round(F, 2), p = signif(p, 3)), row.names = FALSE)
    cat('random-effect SDs: ', paste(sprintf('%s %.4f', vc$grp, vc$sdcor), collapse = ', '), if (lme4::isSingular(f)) '(singular)' else '', '\n')
    cat('slope per shared gene, by input (95% CI):\n')
    print(g_slope[[s_]] |> transmute(input, slope = signif(slope, 3), lo = signif(lo, 3), hi = signif(hi, 3), p = signif(p, 3)), row.names = FALSE)
    cat('slope differences (unadjusted):\n')
    print(g_pair[[s_]] |> transmute(contrast, diff = signif(diff, 3), lo = signif(lo, 3), hi = signif(hi, 3), df = round(df, 1), p = signif(p, 3)), row.names = FALSE)
    em_ends = as.data.frame(emmeans::emmeans(f, ~ input | n, at = list(n = c(genes_ref, round(0.2 * genes_ref))), lmer.df = 'kenward-roger'))
    cat(sprintf('fitted improvement (mean over tissues and algorithms) at %d / %d shared genes: %s\n', genes_ref, round(0.2 * genes_ref),
                paste(sapply(levels(x$input), \(i) sprintf('%s %+.3f / %+.3f', i, em_ends$emmean[em_ends$input == i][1], em_ends$emmean[em_ends$input == i][2])), collapse = '; ')))
    # differences between inputs along the number of shared genes (all pairs, mean over tissues, unadjusted 95% CI)
    g_dline[[s_]] = as.data.frame(summary(pairs(emmeans::emmeans(f, ~ input | n, at = list(n = seq(max(x$n), min(x$n), length.out = 40)),
                                                                 lmer.df = 'kenward-roger'), adjust = 'none'), infer = TRUE)) |>
      transmute(score = s_, n, contrast = gsub('[()]', '', contrast), diff = estimate, lo = lower.CL, hi = upper.CL)
    # inputs compared at low_genes shared genes: all pairs, averaged over tissues and within each tissue (unadjusted)
    lowc = function(em) as.data.frame(summary(pairs(em, adjust = 'none'), infer = TRUE))
    g_low[[s_]] = bind_rows(lowc(emmeans::emmeans(f, ~ input, at = list(n = low_genes), lmer.df = 'kenward-roger')) |> mutate(tissue = 'all'),
                            lowc(emmeans::emmeans(f, ~ input | tissue, at = list(n = low_genes), lmer.df = 'kenward-roger')) |> mutate(tissue = as.character(tissue))) |>
      transmute(score = s_, tissue, contrast = gsub('[()]', '', contrast), diff = estimate, lo = lower.CL, hi = upper.CL, df, p = p.value)
    cat(sprintf('inputs compared at %d shared genes (unadjusted):\n', low_genes))
    print(g_low[[s_]] |> transmute(tissue = ifelse(tissue == 'all', 'mean over tissues', tissue_labels[tissue]), contrast,
                                   diff = sprintf('%+.3f [%+.3f, %+.3f]', diff, lo, hi), df = round(df, 1), p = signif(p, 3)), row.names = FALSE)
    # common-slope check: the same slopes estimated within each tissue
    g_tissue[[s_]] = map_dfr(levels(x$tissue), \(t) {
      ft = suppressMessages(lmerTest::lmer(impr ~ input * n + algorithm * input + (1 | dataset) + (1 | dataset:input), data = filter(x, tissue == t),
                                           contrasts = list(input = 'contr.sum', algorithm = 'contr.sum')))
      as.data.frame(summary(emmeans::emtrends(ft, ~ input, var = 'n', lmer.df = 'kenward-roger'), infer = TRUE)) |>
        transmute(score = s_, tissue = t, input, slope = n.trend, lo = lower.CL, hi = upper.CL, p = p.value)
    })
    cat('slope per shared gene within each tissue (check of the common slope):\n')
    print(g_tissue[[s_]] |> select(tissue, input, slope) |> mutate(slope = signif(slope, 3)) |> pivot_wider(names_from = input, values_from = slope) |>
            as.data.frame(), row.names = FALSE)
  }
  g_anova = bind_rows(g_anova); g_slope = bind_rows(g_slope); g_pair = bind_rows(g_pair); g_tissue = bind_rows(g_tissue); g_fit = bind_rows(g_fit); g_low = bind_rows(g_low); g_dline = bind_rows(g_dline)
  write_csv(bind_rows(g_anova |> mutate(table = 'anova'), g_slope |> mutate(table = 'slope_by_input', input = as.character(input)),
                      g_pair |> mutate(table = 'slope_difference'), g_tissue |> mutate(table = 'slope_by_tissue', input = as.character(input)),
                      g_low |> mutate(table = sprintf('inputs_at_%d_genes', low_genes))) |>
              relocate(table), file.path(out, paste0('input_genes_model', suffix, '.csv')))

  # ---- figure: (A) fitted improvement by input against shared genes, for overall, bio and batch; (B) differences between
  # inputs along the number of shared genes; appendix figure: slopes within each tissue ----
  sc_f = \(v) factor(score_labels[v], levels = score_labels)
  pts = gd |> group_by(score, tissue, n, input) |> summarise(impr = mean(impr), .groups = 'drop') |> mutate(score = sc_f(score))
  fo = g_fit |> mutate(score = sc_f(score))
  ann = g_anova |> filter(term == 'input:n') |>
    transmute(score = sc_f(score), label = sprintf('Input × shared orthologues:\nF(%d, %.0f) = %.2f, p %s', NumDF, DenDF, F,
                                                  ifelse(p < 0.001, '< 0.001', sprintf('= %.3f', p))))
  gA = ggplot(fo, aes(n, emmean, colour = input, fill = input)) +
    geom_hline(yintercept = 0, colour = 'grey40', linetype = 'dashed') +
    geom_point(data = pts, aes(n, impr, colour = input, shape = tissue), inherit.aes = FALSE, size = 1.7, alpha = 0.5) +
    geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.05) +
    facet_wrap(~ score, nrow = 1, scales = 'free_y') +
    scale_colour_manual(values = in_cols, name = 'Input') + scale_fill_manual(values = in_cols, name = 'Input') +
    scale_shape_manual(values = tissue_shapes, labels = tissue_labels, name = 'Tissue (mean over algorithms)') +
    scale_x_reverse(breaks = seq(0, 125, 25)) +
    labs(x = 'Number of shared orthologues', y = 'Improvement over Unintegrated') +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = 'bottom', legend.box = 'horizontal',
                                     strip.background = element_rect(fill = 'grey95', colour = 'grey70'), strip.text = element_text(face = 'bold'))
  # (B) difference between each pair of inputs (first minus second) along the number of shared genes: fitted (mean over tissues,
  # 95% CI) with the tissue x dataset means (points); the slope of each line is the difference in slope between the two inputs,
  # in the printout and CSV (slope_difference; not annotated on the figure). emmeans gives a - b in input order; flipped to b - a.
  lab_d = c(`Native - UINMF` = 'UINMF minus native', `Native - scVI-latent` = 'scVI-latent minus native',
            `scVI-latent - UINMF` = 'UINMF minus scVI-latent')
  fl = \(d) d |> mutate(contrast = factor(lab_d[contrast], levels = lab_d), score = sc_f(score))
  dl = g_dline |> mutate(diff = -diff, new_lo = -hi, new_hi = -lo, lo = new_lo, hi = new_hi) |> fl()
  dp = gd |> group_by(score, tissue, n, input) |> summarise(impr = mean(impr), .groups = 'drop') |> mutate(input = as.character(input)) |>
    pivot_wider(names_from = input, values_from = impr) |>
    transmute(score, tissue = factor(tissue, levels = names(tissue_labels)), n,
              `Native - UINMF` = UINMF - Native, `Native - scVI-latent` = `scVI-latent` - Native, `scVI-latent - UINMF` = UINMF - `scVI-latent`) |>
    pivot_longer(c(`Native - UINMF`, `Native - scVI-latent`, `scVI-latent - UINMF`), names_to = 'contrast', values_to = 'diff') |> fl()
  da = g_pair |> filter(contrast %in% names(lab_d)) |>
    transmute(score, contrast, label = sprintf('Change in the difference per orthologue removed:\n%+.5f (p %s)', diff, ifelse(p < 0.001, '< 0.001', sprintf('= %.3f', p)))) |> fl()
  gB = ggplot(dl, aes(n, diff)) +
    geom_hline(yintercept = 0, colour = 'grey40', linetype = 'dashed') +
    geom_point(data = dp, aes(colour = tissue, shape = tissue), size = 1.7, alpha = 0.6) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = 'grey60', alpha = 0.3) +
    geom_line(linewidth = 1, colour = 'black') +
    facet_grid(contrast ~ score, scales = 'free_y') +
    scale_colour_manual(values = tissue_cols, labels = tissue_labels, name = NULL) +
    scale_shape_manual(values = tissue_shapes, labels = tissue_labels, name = NULL) +
    scale_x_reverse(breaks = seq(0, 125, 25)) +
    labs(x = 'Number of shared orthologues',
         y = 'Difference in improvement (first input minus second)\nline: mean over tissues and algorithms, 95% CI; points: each tissue and dataset') +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = 'bottom',
                                     strip.background = element_rect(fill = 'grey95', colour = 'grey70'), strip.text = element_text(face = 'bold'))
  pG = (gA / gB) + patchwork::plot_layout(heights = c(1, 1.45)) +
    patchwork::plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 13))
  # appendix: slopes within each tissue (check of the common slope)
  gS = ggplot(g_tissue |> mutate(score = sc_f(score), input = factor(input, levels = rev(c('Native', inputs))),
                                 tissue = factor(tissue, levels = names(tissue_labels))), aes(slope, input)) +
    geom_vline(xintercept = 0, colour = 'grey40') +
    geom_point(aes(colour = tissue, shape = tissue), size = 2.4) +
    geom_point(data = g_slope |> mutate(score = sc_f(score), input = factor(input, levels = rev(c('Native', inputs)))),
               aes(slope, input), shape = 124, size = 6, colour = 'black') +
    facet_wrap(~ score, nrow = 1, scales = 'free_x') + scale_x_continuous(labels = \(v) format(v, scientific = FALSE, drop0trailing = TRUE)) +
    scale_colour_manual(values = tissue_cols, labels = tissue_labels, name = NULL) +
    scale_shape_manual(values = tissue_shapes, labels = tissue_labels, name = NULL) +
    labs(x = 'Change in improvement per shared orthologue, within each tissue (bar = common slope from the model)', y = NULL) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), legend.position = 'bottom',
                                     strip.background = element_rect(fill = 'grey95', colour = 'grey70'), strip.text = element_text(face = 'bold'),
                                     plot.margin = margin(5.5, 18, 5.5, 5.5))
  for (ext in c('pdf', 'png')) {
    ggsave(file.path(out, paste0('input_genes_model', suffix, '.', ext)), pG, width = 15, height = 15, dpi = 300, bg = 'white',
           device = if (ext == 'pdf') cairo_pdf else 'png')
    ggsave(file.path(out, paste0('input_genes_model_tissue_slopes', suffix, '.', ext)), gS, width = 13, height = 4.2, dpi = 300, bg = 'white',
           device = if (ext == 'pdf') cairo_pdf else 'png')
  }
}
