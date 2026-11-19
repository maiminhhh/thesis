# Label sensitivity of the 05 method rankings, with the comparison labels from analysis/best_annotator_select.R:
# per dataset and species, the annotator (of the four; no vote, no EBCC) with the highest mean marker AUROC
# on the SELECTION half of cells. Bone marrow, fetal liver and colon were rescored for this
# (LABEL_SET=best_select, scores in label_best_select/); breast cancer's earlier best_annotator rescoring
# used the same columns (human Seurat, mouse Seurat (cross)), so it is reused. Colon human uses SingleR,
# with its pruned calls filled by same-species Seurat (colon/utils.py: 'celltype_singleR|celltype_seurat').
# Per dataset: Kendall's tau between the method rankings under EBCC and under these labels, for the overall
# score (0.6 bio + 0.4 batch) and for bio, and the top 3 under each (as fetal_liver/label_sensitivity.R).
# Usage: Rscript analysis/label_sensitivity_select.R   (from the repository root)
#   -> <tissue>/data/05_evaluated_score/label_sensitivity_select{,_scores}.csv, figures/cross_tissue/label_sensitivity_select.csv
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(readr); library(stringr) })
setwd(Sys.getenv('PROJECT_ROOT', '.'))

batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label',
                  'isolated_label_F1', 'isolated_label_silhouette', 'cLISI')
dataset_names = c('full_mouse', 'dissected_mouse', 'genes2', 'genes4', 'genes6', 'genes8')
tissue_cfg = list(
  colon         = list(label_set = 'best_select',    layout = 'eval'),
  bone_marrow   = list(label_set = 'best_select',    layout = 'eval'),
  fetal_liver   = list(label_set = 'best_select',    layout = 'eval'),
  breast_cancer = list(label_set = 'best_annotator', layout = 'data'))
score_dir = function(tissue, layout, d) {
  if (layout == 'data' || startsWith(d, 'genes')) sprintf('%s/data/%s/', tissue, d)
  else sprintf('%s/data/05_evaluated_score/%s/', tissue, d)
}
best = read_csv('figures/cross_tissue/best_annotator_select.csv', show_col_types = FALSE) |> filter(best)

set_scores = function(dir) {
  files = list.files(dir, pattern = '^05\\..*\\.csv$', full.names = TRUE)
  if (length(files) == 0) return(NULL)
  map_dfr(files, function(f) {
    df = read_csv(f, show_col_types = FALSE)
    tibble(metric = df[[1]], value = as.numeric(df[[2]]), method = names(df)[2])
  }) |>
    mutate(group = case_when(metric %in% batch_metrics ~ 'batch', metric %in% bio_metrics ~ 'bio')) |>
    filter(!is.na(group), !is.na(value)) |>
    group_by(method, group) |> summarise(score = mean(value), .groups = 'drop') |>
    pivot_wider(names_from = group, values_from = score) |>
    mutate(overall = 0.6 * bio + 0.4 * batch)
}
top3 = function(s, col) s |> arrange(desc(.data[[col]])) |> head(3) |> pull(method) |> paste(collapse = ', ')
same = \(x, y) setequal(str_split_1(x, ', '), str_split_1(y, ', '))

all = imap_dfr(tissue_cfg, function(cfg, tissue) {
  alt_desc = best |> filter(tissue == !!tissue) |> arrange(species) |>
    with(paste(sprintf('%s: %s', species, annotator), collapse = '; '))
  scores = map_dfr(dataset_names, function(d) {
    dir = score_dir(tissue, cfg$layout, d)
    e = set_scores(dir); a = set_scores(file.path(dir, str_glue('label_{cfg$label_set}')))
    if (is.null(e) || is.null(a)) { message(tissue, ' ', d, ': no scores under one label set'); return(NULL) }
    bind_rows(mutate(e, labels = 'ebcc'), mutate(a, labels = 'best_select')) |> mutate(dataset = d)
  })
  sens = map_dfr(unique(scores$dataset), function(d) {
    e = filter(scores, dataset == d, labels == 'ebcc'); a = filter(scores, dataset == d, labels == 'best_select')
    m = inner_join(e, a, by = 'method', suffix = c('_ebcc', '_alt'))
    missing = setdiff(union(e$method, a$method), m$method)
    if (length(missing)) message(tissue, ' ', d, ': scored under only one label set: ', paste(missing, collapse = ', '))
    e = filter(e, method %in% m$method); a = filter(a, method %in% m$method)
    tibble(tissue = tissue, alternative = alt_desc, dataset = d, n_methods = nrow(m),
           tau_overall = cor(m$overall_ebcc, m$overall_alt, method = 'kendall'),
           tau_bio     = cor(m$bio_ebcc, m$bio_alt, method = 'kendall'),
           top3_overall_ebcc = top3(e, 'overall'), top3_overall_alt = top3(a, 'overall'),
           top3_bio_ebcc = top3(e, 'bio'), top3_bio_alt = top3(a, 'bio')) |>
      mutate(same_top3_overall = same(top3_overall_ebcc, top3_overall_alt),
             same_top3_bio = same(top3_bio_ebcc, top3_bio_alt))
  })
  write_csv(sens, sprintf('%s/data/05_evaluated_score/label_sensitivity_select.csv', tissue))
  write_csv(scores, sprintf('%s/data/05_evaluated_score/label_sensitivity_select_scores.csv', tissue))
  sens
})
write_csv(all, 'figures/cross_tissue/label_sensitivity_select.csv')
print(as.data.frame(select(all, tissue, dataset, n_methods, tau_overall, tau_bio, same_top3_overall, same_top3_bio)))
message('label sensitivity (selection-half best annotator) done')
