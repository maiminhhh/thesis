# Section 4 of ../EBCC_PIPELINE.md: make the saved fit for the chosen frac (ebcc_frac_choose.R) the
# final fit, without refitting: runs the 02 script's EBCC setup, loads the saved fit as `fit`, runs the
# script's lines after its run_ebcc call (diagnostics, save, dotplots), then writes M and frac back.
# usage: Rscript ebcc_frac_finish.R human|mouse
which = commandArgs(TRUE)[1]
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
pdf(NULL)
chosen = read.csv(sprintf('data/02_annotated_data/%s_ebcc_frac_auroc_chosen.csv', which))
file = if (which == 'human') '02_human_annotated.R' else '02_mouse_annotate.qmd'
src = readLines(file)
find1 = function(pattern) { i = which(trimws(src) == pattern); if (length(i) != 1) stop(file, ': expected one line "', pattern, '"'); i }
setup = find1("setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))")
ebcc  = find1("print('EBCC')")
sweep = find1("Ms = c(1, 2, 3, 5, 10)")
fit_e = which(grepl('^fit = run_ebcc\\(', trimws(src))); stopifnot(length(fit_e) == 1)
fences = which(trimws(src) == '```')
end = if (length(fences)) min(fences[fences > fit_e]) - 1 else length(src)

block = src[ebcc:(sweep - 1)]
i = grep('^frac = ', trimws(block)); stopifnot(length(i) == 1)
block[i] = sprintf('frac = %s', format(chosen$frac))

eval(parse(text = src[setup:(setup + 3)]), envir = globalenv())
eval(parse(text = block), envir = globalenv())
M = chosen$M
fit = readRDS(sprintf('data/02_annotated_data/ebcc_frac_fits/%s_frac%s.rds', which, format(chosen$frac)))
stopifnot(fit$M == M, fit$frac == frac)
message(sprintf('\n===== %s: lines after the run_ebcc call, with the saved fit M = %g, frac = %g [%s] =====',
                file, M, frac, format(Sys.time())))
eval(parse(text = src[(fit_e + 1):end]), envir = globalenv())

# write the choice back into the 02 script
i_M = which(grepl('^M = [0-9]', trimws(src))); i_f = which(grepl('^frac = [0-9]', trimws(src)))
stopifnot(length(i_M) == 1, length(i_f) == 1)
src[i_M] = sprintf('M = %g   # chosen by the held-out (M, frac) grid, 1-SE rule (see %s_ebcc_heldout_grid.csv)', M, which)
src[i_f] = sprintf('frac = %g   # chosen by marker AUROC with M fixed, 50/50 split (see %s_ebcc_frac_auroc_chosen.csv)', frac, which)
writeLines(src, file)
message(sprintf('final EBCC fit saved with M = %g, frac = %g; written back into %s [%s]', M, frac, file, format(Sys.time())))
