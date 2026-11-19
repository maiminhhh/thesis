# Label sensitivity of the 05 method rankings (../EBCC_PIPELINE.md section 9b, targeted check).
# 05 is rescored once with the best-annotator labels (per species the annotator with the highest
# mean marker AUROC, EBCC if it is best; data/05_evaluated_score/label_sensitivity_choice.json):
#   LABEL_SET=best_annotator LABEL_COLS=... bash run_05.sh   -> <scores>/label_best_annotator/
# Per dataset, the method rankings under EBCC and under the best-annotator labels are compared
# with Kendall's tau, for the overall score (0.6 bio + 0.4 batch) and for bio, along with the
# top-3 methods under each. scANVI is trained on the EBCC labels, so it is flagged as label-informed.
#   data/05_evaluated_score/label_sensitivity.csv          tau and top-3 per dataset
#   data/05_evaluated_score/label_sensitivity_scores.csv   every method's scores under both label sets
#   ../figures/fetal_liver/scores/label_sensitivity.png    EBCC vs best-annotator overall score
# usage: Rscript label_sensitivity.R
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'fetal_liver'))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2); library(stringr)
})

label_set = 'best_annotator'
choice = jsonlite::read_json('data/05_evaluated_score/label_sensitivity_choice.json')
alt_desc = paste(map_chr(choice$choice, \(x) str_glue('{x$species}: {x$annotator}')), collapse = ', ')
datasets = c(full_mouse      = 'data/05_evaluated_score/full_mouse/',
             dissected_mouse = 'data/05_evaluated_score/dissected_mouse/',
             genes2 = 'data/genes2/', genes4 = 'data/genes4/',
             genes6 = 'data/genes6/', genes8 = 'data/genes8/')
fig_dir = '../figures/fetal_liver/scores'

# same metric groups as visualisations.R
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label',
                  'isolated_label_F1', 'isolated_label_silhouette', 'cLISI')

read_metrics = function(f) {
  df = read_csv(f, show_col_types = FALSE)
  tibble(metric = df[[1]], value = as.numeric(df[[2]]), method = names(df)[2])
}

# per method: mean of the available batch / bio metrics, overall = 0.6 bio + 0.4 batch
set_scores = function(dir) {
  files = list.files(dir, pattern = '^05\\..*\\.csv$', full.names = TRUE)
  if (length(files) == 0) return(NULL)
  map_dfr(files, read_metrics) |>
    mutate(group = case_when(metric %in% batch_metrics ~ 'batch',
                             metric %in% bio_metrics   ~ 'bio')) |>
    filter(!is.na(group), !is.na(value)) |>
    group_by(method, group) |>
    summarise(score = mean(value), .groups = 'drop') |>
    pivot_wider(names_from = group, values_from = score) |>
    mutate(overall = 0.6 * bio + 0.4 * batch)
}

scores = imap_dfr(datasets, function(dir, d) {
  map_dfr(c('ebcc', label_set), function(ls) {
    s = set_scores(if (ls == 'ebcc') dir else file.path(dir, str_glue('label_{ls}')))
    if (is.null(s)) { message('no scores: ', d, ' ', ls); return(NULL) }
    mutate(s, dataset = d, label_set = ls)
  })
}) |>
  mutate(label_informed = method == 'scANVI')
stopifnot(nrow(scores) > 0)

top3 = function(s, col) s |> arrange(desc(.data[[col]])) |> head(3) |> pull(method) |> paste(collapse = ', ')

sensitivity = map_dfr(names(datasets), function(d) {
  e = filter(scores, dataset == d, label_set == 'ebcc')
  a = filter(scores, dataset == d, label_set == !!label_set)
  if (nrow(e) == 0 || nrow(a) == 0) return(NULL)
  m = inner_join(e, a, by = 'method', suffix = c('_ebcc', '_alt'))
  missing = setdiff(union(e$method, a$method), m$method)
  if (length(missing)) message(d, ': scored under only one label set: ', paste(missing, collapse = ', '))
  e = filter(e, method %in% m$method); a = filter(a, method %in% m$method)
  tibble(dataset = d, label_set = label_set, n_methods = nrow(m),
         tau_overall = cor(m$overall_ebcc, m$overall_alt, method = 'kendall'),
         tau_bio     = cor(m$bio_ebcc, m$bio_alt, method = 'kendall'),
         top3_overall_ebcc = top3(e, 'overall'), top3_overall_alt = top3(a, 'overall'),
         top3_bio_ebcc     = top3(e, 'bio'),     top3_bio_alt      = top3(a, 'bio'))
}) |>
  mutate(same_top3_overall = map2_lgl(top3_overall_ebcc, top3_overall_alt,
                                      \(x, y) setequal(str_split_1(x, ', '), str_split_1(y, ', '))),
         same_top3_bio = map2_lgl(top3_bio_ebcc, top3_bio_alt,
                                  \(x, y) setequal(str_split_1(x, ', '), str_split_1(y, ', '))))

write_csv(sensitivity, 'data/05_evaluated_score/label_sensitivity.csv')
write_csv(scores, 'data/05_evaluated_score/label_sensitivity_scores.csv')
print(as.data.frame(select(sensitivity, dataset, n_methods, tau_overall, tau_bio,
                           same_top3_overall, same_top3_bio)))

# EBCC vs best-annotator overall score, one point per method, tau in each panel
plot_df = scores |>
  filter(label_set != 'ebcc') |>
  inner_join(filter(scores, label_set == 'ebcc') |> select(dataset, method, overall_ebcc = overall),
             by = c('dataset', 'method')) |>
  mutate(dataset = factor(dataset, levels = names(datasets)))
tau_df = sensitivity |>
  mutate(dataset = factor(dataset, levels = names(datasets)),
         text = sprintf('tau overall = %.2f\ntau bio = %.2f', tau_overall, tau_bio))
p = ggplot(plot_df, aes(overall_ebcc, overall)) +
  geom_abline(linetype = 'dashed', colour = 'grey50') +
  geom_point(aes(shape = ifelse(label_informed, 'scANVI (trained on EBCC labels)', 'other methods')),
             colour = 'grey20', size = 2) +
  ggrepel::geom_text_repel(aes(label = method), size = 2.2, max.overlaps = 20, seed = 1) +
  geom_text(data = tau_df, aes(x = -Inf, y = Inf, label = text),
            hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ dataset, nrow = 2) +
  scale_shape_manual(values = c('other methods' = 16, 'scANVI (trained on EBCC labels)' = 17), name = NULL) +
  labs(x = 'Overall score, EBCC labels (0.6 bio + 0.4 batch)',
       y = 'Overall score, best-annotator labels',
       title = 'fetal_liver: method scores under EBCC vs best-annotator labels',
       subtitle = str_glue('best-annotator labels ({alt_desc}): annotator with the highest mean marker AUROC per species')) +
  theme_bw() + theme(legend.position = 'bottom')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(fig_dir, 'label_sensitivity.png'), p, width = 14, height = 9, dpi = 200)
print('label sensitivity done')
