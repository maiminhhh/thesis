# Label sensitivity check, step 4 (../EBCC_PIPELINE.md section 9b): method rankings under EBCC vs under the
# best-annotator labels (data/<dataset>/label_best_annotator/05.*.csv). Scores as in visualisations.R:
# batch / bio = mean of the available metrics; overall = 0.6 bio + 0.4 batch.
# Writes data/05_evaluated_score/label_sensitivity.csv, label_sensitivity_scores.csv and
# ../figures/breast_cancer/scores/label_sensitivity.png. scANVI is trained on EBCC labels (label-informed).
suppressMessages({library(tidyverse); library(jsonlite)})
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'breast_cancer'))

datasets  = c('full_mouse', 'dissected_mouse', 'genes2', 'genes4', 'genes6', 'genes8')   # as in datasets.py
label_set = 'best_annotator'
choice = read_json('data/05_evaluated_score/label_sensitivity_choice.json', simplifyVector = TRUE)
alt_name = paste(sprintf('%s: %s', choice$choice$species, choice$choice$annotator), collapse = ', ')
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label',
                  'isolated_label_F1', 'isolated_label_silhouette', 'cLISI')

read_scores = function(dir) {
  files = list.files(dir, pattern = '^05\\..*\\.csv$', full.names = TRUE)
  if (length(files) == 0) return(NULL)
  map_dfr(files, function(f) {
    df = read_csv(f, show_col_types = FALSE)
    tibble(metric = df[[1]], value = as.numeric(df[[2]]), method = names(df)[2])
  }) |>
    mutate(group = case_when(metric %in% batch_metrics ~ 'batch', metric %in% bio_metrics ~ 'bio')) |>
    filter(!is.na(group), !is.na(value)) |>
    group_by(method, group) |>
    summarise(score = mean(value), .groups = 'drop') |>
    pivot_wider(names_from = group, values_from = score) |>
    mutate(overall = 0.6 * bio + 0.4 * batch)
}

scores = map_dfr(datasets, function(d) {
  bind_rows(read_scores(sprintf('data/%s/', d)) |> mutate(labels = 'ebcc'),
            read_scores(sprintf('data/%s/label_%s/', d, label_set)) |> mutate(labels = label_set)) |>
    mutate(dataset = d)
}) |>
  relocate(dataset, labels, method)

top3 = function(df, col) paste(head(df$method[order(-df[[col]])], 3), collapse = '; ')
summary = map_dfr(datasets, function(d) {
  e = filter(scores, dataset == d, labels == 'ebcc')
  a = filter(scores, dataset == d, labels == label_set)
  m = inner_join(e, a, by = 'method', suffix = c('_ebcc', '_alt'))
  if (nrow(m) < 3) return(tibble(dataset = d, n_methods = nrow(m)))
  tibble(dataset = d, alternative = alt_name, n_methods = nrow(m),
         tau_overall = cor(m$overall_ebcc, m$overall_alt, method = 'kendall', use = 'complete.obs'),
         tau_bio     = cor(m$bio_ebcc, m$bio_alt, method = 'kendall', use = 'complete.obs'),
         top3_overall_ebcc = top3(e, 'overall'), top3_overall_alt = top3(a, 'overall'),
         top3_bio_ebcc = top3(e, 'bio'), top3_bio_alt = top3(a, 'bio'),
         same_top3_overall = setequal(head(e$method[order(-e$overall)], 3), head(a$method[order(-a$overall)], 3)),
         same_top3_bio = setequal(head(e$method[order(-e$bio)], 3), head(a$method[order(-a$bio)], 3)))
})
summary$note = 'scANVI is trained on EBCC labels (label-informed) in both scorings'

dir.create('data/05_evaluated_score', showWarnings = FALSE)
write_csv(summary, 'data/05_evaluated_score/label_sensitivity.csv')
write_csv(scores, 'data/05_evaluated_score/label_sensitivity_scores.csv')
print(as.data.frame(summary |> select(dataset, n_methods, tau_overall, tau_bio, same_top3_overall, same_top3_bio)))

plot_df = scores |>
  select(dataset, labels, method, overall) |>
  pivot_wider(names_from = labels, values_from = overall) |>
  filter(!is.na(ebcc), !is.na(.data[[label_set]])) |>
  mutate(dataset = factor(dataset, levels = datasets), uinmf = if_else(grepl('UINMF', method), 'UINMF', 'Base'))
taus = summary |> filter(!is.na(tau_overall)) |>
  mutate(dataset = factor(dataset, levels = datasets), lab = sprintf('tau = %.2f', tau_overall))
p = ggplot(plot_df, aes(ebcc, .data[[label_set]])) +
  geom_abline(linetype = 'dashed', colour = 'grey60') +
  geom_point(aes(shape = uinmf), size = 2) +
  geom_text(data = filter(plot_df, method == 'scANVI'), aes(label = 'scANVI (label-informed)'),
            size = 2.5, hjust = -0.1, colour = 'grey30') +
  geom_text(data = taus, aes(x = -Inf, y = Inf, label = lab), hjust = -0.1, vjust = 1.5, size = 3.2) +
  facet_wrap(~ dataset, nrow = 1) +
  labs(x = 'overall score, EBCC labels (0.6 bio + 0.4 batch)',
       y = 'overall score, best-annotator labels', shape = NULL,
       title = 'breast_cancer: method scores, EBCC vs best-annotator labels',
       subtitle = sprintf('best annotator by marker AUROC: %s', alt_name)) +
  theme_bw()
dir.create('../figures/breast_cancer/scores', showWarnings = FALSE, recursive = TRUE)
ggsave('../figures/breast_cancer/scores/label_sensitivity.png', p, width = 16, height = 4, bg = 'white')
cat('label sensitivity done\n')
