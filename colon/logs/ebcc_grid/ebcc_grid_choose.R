# Combine the (M, frac) held-out grid for one species, pick (M, frac) and write it back into the 02 script.
# Rule: among combinations whose mean held-out ll is within 1 SE (SE of the best) of the best mean,
# take the smallest M, then the smallest frac.
# usage: Rscript ebcc_grid_choose.R human|mouse      (prints "CHANGED yes|no" as its last line)
# ("changed" compares with the M and frac in the script; the labels themselves changed with the
#  reference fix of 2026-10-23, so the final fit is rerun either way)
which = commandArgs(TRUE)[1]
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
suppressMessages(library(ggplot2))
file = if (which == 'human') '02_human_annotated.R' else '02_mouse_annotate.qmd'
src = readLines(file)

parts = lapply(list.files('logs/ebcc_grid', sprintf('^%s_frac.*\\.rds$', which), full.names = TRUE), readRDS)
fracs = sapply(parts, `[[`, 'frac')
stopifnot(!anyDuplicated(fracs), setequal(fracs, c(0, 0.25, 0.5, 1, 2)),
          all(sapply(parts, function(p) identical(p$trusted, parts[[1]]$trusted))),
          all(sapply(parts, function(p) identical(p$ho, parts[[1]]$ho))))
Ms = c(1, 2, 3, 5, 10)

grid = do.call(rbind, lapply(parts, function(p) data.frame(
  M = Ms, frac = p$frac, mean_ll = colMeans(p$ll),
  se = apply(p$ll, 2, sd) / sqrt(nrow(p$ll)), n_heldout = nrow(p$ll))))
grid = grid[order(grid$M, grid$frac), ]
best = which.max(grid$mean_ll)
grid$within_1se_of_best = grid$mean_ll >= grid$mean_ll[best] - grid$se[best]
grid$best = seq_len(nrow(grid)) == best
cand = grid[grid$within_1se_of_best, ]
pick = cand[order(cand$M, cand$frac), ][1, ]
grid$chosen = grid$M == pick$M & grid$frac == pick$frac
rownames(grid) = NULL
print(grid, digits = 6)

# previous choice, as currently written in the script
prev_M    = as.numeric(sub('^M = ([0-9.]+).*', '\\1', grep('^M = [0-9]', trimws(src), value = TRUE)))
prev_frac = as.numeric(sub('^frac = ([0-9.]+).*', '\\1', grep('^frac = [0-9]', trimws(src), value = TRUE)))
stopifnot(length(prev_M) == 1, length(prev_frac) == 1)
changed = pick$M != prev_M || pick$frac != prev_frac

rule = paste('mean +- SE (sd/sqrt(n)) of the per-cell held-out log-likelihood for each (M, frac);',
             'among combinations within 1 SE of the best mean, the smallest M, then the smallest frac')
out = sprintf('data/02_annotated_data/%s_ebcc_heldout_grid', which)
saveRDS(list(ll = setNames(lapply(parts, `[[`, 'll'), fracs), summary = grid, rule = rule,
             chosen = c(M = pick$M, frac = pick$frac), trusted = parts[[1]]$trusted), paste0(out, '.rds'))
write.csv(grid, paste0(out, '.csv'), row.names = FALSE)
write.csv(data.frame(species = which, M = pick$M, frac = pick$frac, changed = changed,
                     previous_M = prev_M, previous_frac = prev_frac),
          paste0(out, '_chosen.csv'), row.names = FALSE)

p = ggplot(grid, aes(factor(frac), factor(M), fill = mean_ll)) +
  geom_tile() +
  geom_tile(data = grid[grid$within_1se_of_best, ], fill = NA, colour = 'black', linewidth = 0.4) +
  geom_tile(data = grid[grid$chosen, ], fill = NA, colour = 'black', linewidth = 1.5) +
  geom_text(aes(label = sprintf('%.4f\n(%.4f)', mean_ll, se)), size = 3) +
  scale_fill_viridis_c() +
  labs(x = 'frac', y = 'M', fill = 'mean\nheld-out ll',
       title = sprintf('%s EBCC held-out log-likelihood', which),
       subtitle = sprintf('thin outline: within 1 SE of best; thick: chosen (M = %g, frac = %g)', pick$M, pick$frac)) +
  theme_minimal()
dir.create('../figures/colon/annotation', showWarnings = FALSE, recursive = TRUE)
ggsave(sprintf('../figures/colon/annotation/ebcc_heldout_grid_%s.png', which), p, width = 7, height = 5, dpi = 150)

# write the choice back into the 02 script
note = sprintf('   # chosen by the held-out (M, frac) grid, 1-SE rule (see %s_ebcc_heldout_grid.csv)', which)
i_M = which(grepl('^M = [0-9]', trimws(src))); i_f = which(grepl('^frac = [0-9]', trimws(src)))
src[i_M] = sprintf('M = %g%s', pick$M, note)
src[i_f] = sprintf('frac = %g%s', pick$frac, note)
writeLines(src, file)

cat(sprintf('%s: chosen M = %g, frac = %g (previously M = %g, frac = %g)\nCHANGED %s\n',
            which, pick$M, pick$frac, prev_M, prev_frac, if (changed) 'yes' else 'no'))
