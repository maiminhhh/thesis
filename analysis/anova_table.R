# LaTeX ANOVA table for section 2.2 (the block model behind section22<suffix>.pdf, analysis/mcb_analysis.R section 8 / mcb()):
# lm(overall improvement ~ method + tissue), dissected, 19 embedding methods (graph-based excluded, as the section 2.2
# figure) x 4 tissues, one value per cell, so the residual is the
# method x tissue interaction. Balanced design: sequential SS = type III (checked below).
# Usage: R_LIBS_SITE= R_LIBS_USER= <r45>/bin/Rscript --vanilla analysis/anova_table.R -> figures/mcb/anova_table<suffix>.tex
setwd(Sys.getenv('PROJECT_ROOT', '.'))
suppressPackageStartupMessages({ library(knitr); library(readr) })

# ---- config ----
suffix = '_MAP_SAS'                                   # score definition, as analysis/mcb_analysis.R ('' = reference)
in_file = sprintf('figures/benchmark_summary/benchmark_improvement%s.csv', suffix)
exclude_graph_based = TRUE                            # as section22<suffix>.pdf; FALSE = all 23 methods (analysis/mcb_analysis.R section 1)
graph_only = c('BBKNN', 'BBKNN UINMF', 'BBKNN scVI-latent', 'SAMap')   # as analysis/mcb_analysis.R
out_file = sprintf('figures/mcb/anova_table%s%s.tex', suffix, if (exclude_graph_based) '_no_graph' else '')

d = read_csv(in_file, show_col_types = FALSE)
if (exclude_graph_based) d = d[!d$method %in% graph_only, ]
d = data.frame(improvement = d$imp_overall, method = factor(d$method), tissue = factor(d$tissue))
stopifnot(nrow(d) == nlevels(d$method) * 4, nlevels(d$method) == if (exclude_graph_based) 19 else 23, all(table(d$method, d$tissue) == 1))

fit = lm(improvement ~ method + tissue, data = d)
a = anova(fit)
a_rev = anova(lm(improvement ~ tissue + method, data = d))
stopifnot(abs(a['method', 'Sum Sq'] - a_rev['method', 'Sum Sq']) < 1e-12)   # balanced: order does not matter

tab = data.frame(Source = c('Method', 'Tissue', 'Method $\\times$ tissue (residual)'),
                 df = a$Df, SS = a$`Sum Sq`, MS = a$`Mean Sq`, F = a$`F value`, p = a$`Pr(>F)`)
tab$Share = 100 * tab$SS / sum(tab$SS)
tab = rbind(tab, data.frame(Source = 'Total', df = sum(tab$df), SS = sum(tab$SS), MS = NA, F = NA, p = NA, Share = 100))

fmt = function(x, digits) ifelse(is.na(x), '', formatC(x, format = 'f', digits = digits))
tab_out = data.frame(
  Source = tab$Source,
  df = tab$df,
  `Sum of squares` = fmt(tab$SS, 4),
  `Mean square` = fmt(tab$MS, 4),
  `$F$` = fmt(tab$F, 2),
  `$p$` = ifelse(is.na(tab$p), '', ifelse(tab$p < 0.001, '$<$ 0.001', fmt(tab$p, 3))),
  `Share (\\%)` = fmt(tab$Share, 0),
  check.names = FALSE)

latex = kable(tab_out, format = 'latex', booktabs = TRUE, escape = FALSE, align = 'lrrrrrr', linesep = '')
writeLines(latex, out_file)

ss = a$`Sum Sq`
cat(latex, sep = '\n')
cat(sprintf('\nPartial eta^2 (method) = %.2f; residual SD = %.4f; wrote %s\n', ss[1] / (ss[1] + ss[3]), sqrt(a['Residuals', 'Mean Sq']), out_file))
