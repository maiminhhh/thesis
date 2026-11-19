# Appendix table for section 2.2: Hsu's MCB on the overall score and on its two components separately, 19 embedding
# methods (graph-based excluded, as section22<suffix>.pdf) x 4 tissues, outcome = improvement over unintegrated (dissected).
# Uses mcb() from analysis/mcb_analysis.R unchanged (lm(y ~ method + tissue); one-sided Dunnett d; interval [min(0, D - d SE),
# max(0, D + d SE)]; could be best = upper > 0).
# Usage: R_LIBS_SITE= R_LIBS_USER= <r45>/bin/Rscript --vanilla analysis/mcb_component_table.R -> figures/mcb/mcb_components<suffix>.{tex,csv}
setwd(Sys.getenv('PROJECT_ROOT', '.'))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })

# ---- config ----
suffix = '_MAP_SAS'
graph_only = c('BBKNN', 'BBKNN UINMF', 'BBKNN scVI-latent', 'SAMap')   # as analysis/mcb_analysis.R
alpha = 0.05                                                           # as analysis/mcb_analysis.R
out = sprintf('figures/mcb/mcb_components%s', suffix)

# mcb() exactly as in analysis/mcb_analysis.R
ex = parse('analysis/mcb_analysis.R')
eval(ex[[which(sapply(ex, \(e) startsWith(deparse(e)[1], 'mcb = function')))]])

d = read_csv(sprintf('figures/benchmark_summary/benchmark_improvement%s.csv', suffix), show_col_types = FALSE) |>
  filter(!method %in% .env$graph_only) |> select(method, tissue, overall = imp_overall, bio = imp_bio, batch = imp_batch)
stopifnot(n_distinct(d$method) == 19, nrow(d) == 76)
res = lapply(c('overall', 'bio', 'batch'), \(o) { r = mcb(d, o)
  cat(sprintf('%-8s d = %.3f, SE = %.4f, d*SE = %.4f, could be best: %d of %d\n', o, r$d, r$se, r$d * r$se, sum(r$res$could_be_best), nrow(r$res)))
  r$res |> transmute(method, score = o, mean, D, mcb_lo, mcb_hi, could_be_best) })
res = bind_rows(res)
write_csv(res, paste0(out, '.csv'))

ord = res |> filter(score == 'overall') |> arrange(-mean) |> pull(method)
lab = function(m) str_replace(m, '^(CCA|RPCA)', '\\1 (Seurat)') |> str_replace('^scANVI$', 'scANVI$^\\\\dagger$')
cell = function(r) {   # mean, then MCB interval; bold when the method could be the best
  s = sprintf('%.3f & [%.3f, %.3f]', r$mean, r$mcb_lo, r$mcb_hi)
  if (r$could_be_best) sprintf('\\textbf{%.3f} & \\textbf{[%.3f, %.3f]}', r$mean, r$mcb_lo, r$mcb_hi) else s
}
rows = sapply(ord, function(m) paste(lab(m), '&', paste(sapply(c('overall', 'bio', 'batch'),
  \(o) cell(res[res$method == m & res$score == o, ])), collapse = ' & '), '\\\\'))
n_best = sapply(c('overall', 'bio', 'batch'), \(o) sum(res$could_be_best[res$score == o]))
tex = c('\\begin{tabular}{lrrrrrr}', '\\toprule',
        '& \\multicolumn{2}{c}{Overall} & \\multicolumn{2}{c}{Bio-conservation} & \\multicolumn{2}{c}{Batch correction} \\\\',
        '\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}',
        'Method & Mean & MCB interval & Mean & MCB interval & Mean & MCB interval \\\\', '\\midrule',
        gsub('_', '\\\\_', rows), '\\midrule',
        sprintf('Could be best & \\multicolumn{2}{c}{%d of 19} & \\multicolumn{2}{c}{%d of 19} & \\multicolumn{2}{c}{%d of 19} \\\\',
                n_best[1], n_best[2], n_best[3]),
        '\\bottomrule', '\\end{tabular}')
writeLines(tex, paste0(out, '.tex'))
cat('\nCould be best:\n'); for (o in c('overall', 'bio', 'batch'))
  cat(sprintf('  %-8s %s\n', o, paste(res$method[res$score == o & res$could_be_best], collapse = ', ')))
cat('Not best:\n'); for (o in c('overall', 'bio', 'batch'))
  cat(sprintf('  %-8s %s\n', o, paste(res$method[res$score == o & !res$could_be_best], collapse = ', ')))
