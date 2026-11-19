# Section 3b of ../EBCC_PIPELINE.md: one full EBCC fit (R = 10, all cells) for one species and frac,
# with M fixed at the held-out grid's choice. Runs the 02 script's own EBCC setup and fit lines
# (frac and M overridden) and saves the fit, so the chosen one can become the final fit unchanged.
# usage: Rscript ebcc_frac_fit.R human|mouse <frac> [n_cores]
#   frac = 0 with "reuse": assemble the fit from the final M/frac=0 fit of 2026-10-23
#   (ebcc_grid_finish.R; params, posterior and the restart ELBOs from its log) instead of refitting
args = commandArgs(TRUE); which = args[1]; frac_new = as.numeric(args[2])
n_cores = if (length(args) >= 3 && args[3] != 'reuse') as.integer(args[3]) else 5L
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
pdf(NULL)
out = sprintf('data/02_annotated_data/ebcc_frac_fits/%s_frac%s.rds', which, format(frac_new))
M_new = read.csv(sprintf('data/02_annotated_data/%s_ebcc_heldout_grid_chosen.csv', which))$M

if (length(args) >= 3 && args[3] == 'reuse') {
  stopifnot(frac_new == 0)
  log = readLines(sprintf('logs/ebcc_grid/finish_%s.log', which))
  stopifnot(any(log == sprintf('EBCC done with M = %g, frac = 0', M_new) | startsWith(log, sprintf('EBCC done with M = %g, frac = 0 [', M_new))))
  fit_at = grep('lines from "# ---- fit ----"', log)
  r = regmatches(log[-(1:fit_at)], regexec('^restart +([0-9]+) \\| iterations +([0-9]+) \\| ELBO (-?[0-9.]+)', log[-(1:fit_at)]))
  r = do.call(rbind, r[lengths(r) == 4])
  r = r[!duplicated(r[, 2]), , drop = FALSE]
  stopifnot(nrow(r) == 10)
  r = r[order(as.integer(r[, 2])), ]
  post = readRDS(sprintf('data/02_annotated_data/%s_ebcc_posterior.rds', which))
  elbos = as.numeric(r[, 4])
  saveRDS(list(label = factor(colnames(post)[max.col(post, ties.method = 'first')], levels = colnames(post)),
               posterior = post, par = readRDS(sprintf('data/02_annotated_data/%s_ebcc_params.rds', which)),
               elbo = max(elbos), elbo_restarts = elbos, iterations = as.integer(r[which.max(elbos), 3]),
               M = M_new, frac = 0, reused = TRUE), out)
  message('frac = 0 fit assembled from the final fit of ebcc_grid_finish.R -> ', out)
  quit(save = 'no')
}

file = if (which == 'human') '02_human_annotated.R' else '02_mouse_annotate.qmd'
src = readLines(file)
find1 = function(pattern) { i = which(trimws(src) == pattern); if (length(i) != 1) stop(file, ': expected one line "', pattern, '"'); i }
setup = find1("setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))")
ebcc  = find1("print('EBCC')")
sweep = find1("Ms = c(1, 2, 3, 5, 10)")
fit_s = find1("# ---- fit ----")
fit_e = which(grepl('^fit = run_ebcc\\(', trimws(src)))
stopifnot(length(fit_e) == 1, fit_e > fit_s)

block = src[ebcc:(sweep - 1)]
i = grep('^frac = ', trimws(block)); stopifnot(length(i) == 1)
block[i] = sprintf('frac = %s', format(frac_new))
fit_block = src[fit_s:fit_e]
i = grep('^M = ', trimws(fit_block)); stopifnot(length(i) == 1)
fit_block[i] = sprintf('M = %g', M_new)

eval(parse(text = src[setup:(setup + 3)]), envir = globalenv())   # setwd, python, library.R, source.R
# restarts in parallel; run_ebcc seeds each restart so this gives the same fit as n_cores = 1
run_ebcc_serial = run_ebcc
run_ebcc = function(...) run_ebcc_serial(..., n_cores = n_cores)
eval(parse(text = block), envir = globalenv())
message(sprintf('%s: full fit M = %g, frac = %s, R = 10 [%s]', which, M_new, format(frac_new), format(Sys.time())))
eval(parse(text = fit_block), envir = globalenv())
stopifnot(frac == frac_new, M == M_new)
saveRDS(list(label = fit$label, posterior = fit$posterior, par = fit$par, elbo = fit$elbo,
             elbo_restarts = fit$elbo_restarts, iterations = fit$iterations,
             M = M, frac = frac, trusted = trusted, reused = FALSE), out)
message(sprintf('saved %s [%s]', out, format(Sys.time())))
