# Section 9b of ../EBCC_PIPELINE.md, steps 1-2: per species, the annotator with the highest mean marker
# AUROC over the cell types it calls (the heatmap's "Mean (all called types)" row), from
# ../figures/colon/annotation/marker_auroc.csv. EBCC is kept for a species where it has the highest mean
# (ties -> EBCC); otherwise that annotator's 03 obs column is used. This only picks the comparison label
# set for the robustness check: EBCC stays the primary label set.
# Writes data/05_evaluated_score/label_sensitivity_choice.json; prints "LABEL_COLS=..." or "SKIP".
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
suppressMessages(source('../library.R'))
obs_col = c(seurat = 'celltype_seurat', singleR = 'celltype_singleR', seurat_cross = 'celltype_seurat_cross',
            singleR_cross = 'celltype_singleR_cross', ebcc = 'ebcc_label')
species_id = c(human = 'hsapiens', mouse = 'mmusculus')

res = read.csv('../figures/colon/annotation/marker_auroc.csv')
stopifnot(setequal(unique(res$annotator), names(obs_col)), setequal(unique(res$species), names(species_id)))
means = aggregate(auroc ~ species + annotator, res, mean)   # NA rows (type not called / no markers) dropped
choice = lapply(setNames(names(species_id), names(species_id)), function(sp) {
  m = means[means$species == sp, ]
  m = m[order(-m$auroc, m$annotator != 'ebcc'), ]          # ties -> EBCC
  list(mean_auroc = setNames(as.list(m$auroc), m$annotator), best = m$annotator[1],
       use = m$annotator[1], obs_column = unname(obs_col[m$annotator[1]]))
})
print(means[order(means$species, -means$auroc), ], digits = 3)
skip = all(sapply(choice, `[[`, 'use') == 'ebcc')
label_cols = paste(sprintf('%s=%s', species_id, sapply(choice[names(species_id)], `[[`, 'obs_column')), collapse = ',')
jsonlite::write_json(list(tissue = 'colon', source = 'figures/colon/annotation/marker_auroc.csv',
                          rule = 'per species: annotator with the highest mean marker AUROC over the types it calls; EBCC if it is highest (ties -> EBCC)',
                          species = choice, skip = skip, label_set = 'best_annotator',
                          label_cols = if (skip) NULL else label_cols),
                     'data/05_evaluated_score/label_sensitivity_choice.json', auto_unbox = TRUE, pretty = TRUE, digits = NA)
for (sp in names(choice)) cat(sprintf('%s: %s (mean AUROC %.3f; EBCC %.3f)\n', sp, choice[[sp]]$use,
                                      choice[[sp]]$mean_auroc[[choice[[sp]]$use]], choice[[sp]]$mean_auroc[['ebcc']]))
cat(if (skip) 'SKIP\n' else sprintf('LABEL_COLS=%s\n', label_cols))
