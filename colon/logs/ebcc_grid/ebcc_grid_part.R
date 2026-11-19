# One column of the (M, frac) held-out grid: run a 02 script's EBCC setup with frac overridden,
# then the held-out log-likelihood for every M in Ms.
# usage: Rscript ebcc_grid_part.R human|mouse <frac>
args = commandArgs(TRUE); which = args[1]; frac_new = as.numeric(args[2])
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
pdf(NULL)   # the mouse EBCC block draws dot plots; don't let parallel jobs share Rplots.pdf

file = if (which == 'human') '02_human_annotated.R' else '02_mouse_annotate.qmd'
src = readLines(file)
find1 = function(pattern) { i = which(trimws(src) == pattern); if (length(i) != 1) stop(file, ': expected one line "', pattern, '"'); i }
setup = find1("setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))")
ebcc  = find1("print('EBCC')")
sweep = find1("Ms = c(1, 2, 3, 5, 10)")
block = src[ebcc:(sweep - 1)]
f_line = grep('^frac = ', trimws(block))
if (length(f_line) != 1) stop('expected one "frac = " line in the EBCC block')
block[f_line] = sprintf('frac = %s', format(frac_new))

eval(parse(text = src[setup:(setup + 3)]), envir = globalenv())   # setwd, python, library.R, source.R
# restarts in parallel; run_ebcc seeds each restart so this gives the same fits as n_cores = 1
run_ebcc_serial = run_ebcc
run_ebcc = function(...) run_ebcc_serial(..., n_cores = 3L)
eval(parse(text = block), envir = globalenv())
stopifnot(frac == frac_new)

Ms = c(1, 2, 3, 5, 10)
message(sprintf('%s frac = %s: fitting M = %s [%s]', which, format(frac), paste(Ms, collapse = ', '), format(Sys.time())))
ll = sapply(Ms, heldout_ll)
colnames(ll) = Ms
saveRDS(list(ll = ll, frac = frac, Ms = Ms, ho = ho, trusted = trusted),
        sprintf('logs/ebcc_grid/%s_frac%s.rds', which, format(frac)))
message('done [', format(Sys.time()), ']')
