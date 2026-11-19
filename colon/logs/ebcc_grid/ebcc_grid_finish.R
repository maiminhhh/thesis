# Final EBCC fit for one species with the (M, frac) that ebcc_grid_choose.R wrote into the 02 script:
# runs the script's EBCC setup (without the held-out sweep), then its fit/diagnostics/save part.
# usage: Rscript ebcc_grid_finish.R human|mouse
which = commandArgs(TRUE)[1]
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
pdf(NULL)
file = if (which == 'human') '02_human_annotated.R' else '02_mouse_annotate.qmd'
src = readLines(file)
find1 = function(pattern) { i = which(trimws(src) == pattern); if (length(i) != 1) stop(file, ': expected one line "', pattern, '"'); i }
setup = find1("setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))")
ebcc  = find1("print('EBCC')")
sweep = find1("Ms = c(1, 2, 3, 5, 10)")
fit   = find1("# ---- fit ----")
fences = which(trimws(src) == '```')
end = if (length(fences)) min(fences[fences > fit]) - 1 else length(src)

eval(parse(text = src[setup:(setup + 3)]), envir = globalenv())
run_ebcc_serial = run_ebcc
run_ebcc = function(...) run_ebcc_serial(..., n_cores = 10L)   # same fits as n_cores = 1
for (b in list(src[ebcc:(sweep - 1)], src[fit:end])) {
  message(sprintf('\n===== %s: lines from "%s" [%s] =====', file, trimws(b[1]), format(Sys.time())))
  eval(parse(text = b), envir = globalenv())
}
message(sprintf('EBCC done with M = %g, frac = %g [%s]', M, frac, format(Sys.time())))
