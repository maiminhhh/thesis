# Label sensitivity check, step 1-2 (../EBCC_PIPELINE.md section 9b): which labels to compare EBCC against.
# From ../figures/breast_cancer/annotation/marker_auroc.csv (written by ../analysis/visualisations.R), each annotator's
# mean AUROC over the cell types it calls; per species keep EBCC if it is best (ties -> EBCC), otherwise
# use the best annotator's 03 column. Writes data/05_evaluated_score/label_sensitivity_choice.json.
# The AUROC is used only to pick a deliberately tough comparison, not to choose the primary labels.
suppressMessages({library(tidyverse); library(jsonlite)})
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'breast_cancer'))

columns = c(seurat = 'celltype_seurat', singleR = 'celltype_singleR', seurat_cross = 'celltype_seurat_cross',
            singleR_cross = 'celltype_singleR_cross', ebcc = 'ebcc_label')   # as `annotators` in visualisations.R
means = read_csv('../figures/breast_cancer/annotation/marker_auroc.csv', show_col_types = FALSE) |>
  filter(!is.na(auroc)) |>
  group_by(species, annotator) |>
  summarise(mean_auroc = mean(auroc), n_types = n(), .groups = 'drop')

choice = means |>
  group_by(species) |>
  arrange(desc(mean_auroc), annotator != 'ebcc', .by_group = TRUE) |>   # ties -> EBCC
  slice(1) |>
  ungroup() |>
  mutate(column = columns[annotator])
species_order = c('human', 'mouse'); obs_species = c(human = 'hsapiens', mouse = 'mmusculus')
choice = choice[match(species_order, choice$species), ]
skip = all(choice$annotator == 'ebcc')
label_cols = paste0(obs_species[choice$species], '=', choice$column, collapse = ',')

dir.create('data/05_evaluated_score', showWarnings = FALSE)
write_json(list(rule = paste('per species, annotator with the highest mean marker AUROC over the types it calls',
                             '(ties -> EBCC); EBCC kept where it is best; check skipped if EBCC is best in both'),
                mean_auroc = means, choice = select(choice, species, annotator, column, mean_auroc),
                skip = skip, label_set = 'best_annotator', label_cols = label_cols),
           'data/05_evaluated_score/label_sensitivity_choice.json', auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(as.data.frame(means |> arrange(species, desc(mean_auroc))))
cat(sprintf('%s: %s (mean AUROC %.3f)\n', choice$species, choice$annotator, choice$mean_auroc), sep = '')
cat(if (skip) 'SKIP: EBCC is the best annotator in both species\n' else sprintf('LABEL_COLS=%s\n', label_cols))
