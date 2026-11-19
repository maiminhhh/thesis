# Label sensitivity of the 05 method rankings (section 9b of ../EBCC_PIPELINE.md).
# The 05 scores under EBCC (primary labels) are compared with the rescoring under the best-annotator
# label set (label_choice.R; LABEL_SET=best_annotator LABEL_COLS=... bash run_05.sh, label_<set>/ subfolders):
# per dataset, Kendall tau between the methods' scores under EBCC and under the alternative, for the
# overall score (0.6 bio + 0.4 batch) and for bio, plus the top-3 methods under each.
# Batch / bio = mean of the available metrics, as in ../analysis/visualisations.R.
# scANVI is trained on ebcc_label (label-informed), so it stays EBCC-informed under every label set.
# Writes data/05_evaluated_score/label_sensitivity{,_scores}.csv, ../figures/colon/scores/label_sensitivity.png
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
suppressMessages(source('../library.R'))

datasets = c(full_mouse = 'data/05_evaluated_score/full_mouse/', dissected_mouse = 'data/05_evaluated_score/dissected_mouse/',
             genes2 = 'data/genes2/', genes4 = 'data/genes4/', genes6 = 'data/genes6/', genes8 = 'data/genes8/')
label_sets = c(ebcc = '', best_annotator = 'label_best_annotator/')
choice = jsonlite::read_json('data/05_evaluated_score/label_sensitivity_choice.json')
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn')
bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label',
                  'isolated_label_F1', 'isolated_label_silhouette', 'cLISI')

read_scores = function(dir) {
  files = list.files(dir, pattern = '^05\\..*\\.csv$', full.names = TRUE)
  if (!length(files)) stop('no 05 scores in ', dir)
  do.call(rbind, lapply(files, function(f) {
    df = read.csv(f, check.names = FALSE)
    v = setNames(as.numeric(df[[2]]), df[[1]])
    data.frame(method = names(df)[2], batch = mean(v[batch_metrics], na.rm = TRUE),
               bio = mean(v[bio_metrics], na.rm = TRUE))
  }))
}

scores = do.call(rbind, lapply(names(datasets), function(d) do.call(rbind, lapply(names(label_sets), function(l) {
  s = read_scores(paste0(datasets[[d]], label_sets[[l]]))
  s$overall = 0.6 * s$bio + 0.4 * s$batch
  s$rank_overall = rank(-s$overall, ties.method = 'min')
  s$rank_bio = rank(-s$bio, ties.method = 'min')
  cbind(dataset = d, label_set = l, s)
}))))
scores$label_informed = scores$method == 'scANVI'
if (!any(scores$label_informed)) warning("no scANVI scores found")
write.csv(scores, 'data/05_evaluated_score/label_sensitivity_scores.csv', row.names = FALSE)

top3 = function(s, col) paste(s$method[order(-s[[col]])][1:3], collapse = '; ')
sens = do.call(rbind, lapply(names(datasets), function(d) do.call(rbind, lapply(setdiff(names(label_sets), 'ebcc'), function(l) {
  e = scores[scores$dataset == d & scores$label_set == 'ebcc', ]
  a = scores[scores$dataset == d & scores$label_set == l, ]
  m = intersect(e$method, a$method)
  if (length(m) < nrow(e) || length(m) < nrow(a))
    warning(d, ' / ', l, ': methods missing from one label set: ', paste(setdiff(union(e$method, a$method), m), collapse = ', '))
  e = e[match(m, e$method), ]; a = a[match(m, a$method), ]
  data.frame(dataset = d, alternative = l, n_methods = length(m),
             tau_overall = cor(e$overall, a$overall, method = 'kendall'),
             tau_bio = cor(e$bio, a$bio, method = 'kendall'),
             top3_ebcc_overall = top3(e, 'overall'), top3_alt_overall = top3(a, 'overall'),
             same_top3_overall = setequal(e$method[order(-e$overall)][1:3], a$method[order(-a$overall)][1:3]),
             top3_ebcc_bio = top3(e, 'bio'), top3_alt_bio = top3(a, 'bio'),
             same_top3_bio = setequal(e$method[order(-e$bio)][1:3], a$method[order(-a$bio)][1:3]))
}))))
write.csv(sens, 'data/05_evaluated_score/label_sensitivity.csv', row.names = FALSE)
print(sens[c('dataset', 'alternative', 'n_methods', 'tau_overall', 'tau_bio', 'same_top3_overall', 'same_top3_bio')], digits = 3)

# EBCC vs alternative overall score, one point per method, tau in each panel
wide = merge(scores[scores$label_set == 'ebcc', c('dataset', 'method', 'overall', 'label_informed')],
             scores[scores$label_set != 'ebcc', c('dataset', 'label_set', 'method', 'overall')],
             by = c('dataset', 'method'), suffixes = c('_ebcc', '_alt'))
names(wide)[names(wide) == 'label_set'] = 'alternative'
wide$dataset = factor(wide$dataset, levels = names(datasets))
sens$dataset = factor(sens$dataset, levels = names(datasets))
sens$tau_label = sprintf('tau = %.2f (bio %.2f)', sens$tau_overall, sens$tau_bio)
p = ggplot(wide, aes(overall_ebcc, overall_alt)) +
  geom_abline(linetype = 'dashed', colour = 'grey60') +
  geom_point(aes(shape = label_informed), size = 1.8) +
  ggrepel::geom_text_repel(aes(label = method), size = 2, max.overlaps = 20, segment.size = 0.2) +
  geom_text(data = sens, aes(x = -Inf, y = Inf, label = tau_label), hjust = -0.05, vjust = 1.3, size = 3) +
  facet_wrap(~ dataset, nrow = 2) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c('no', 'yes (scANVI, trained on EBCC)')) +
  labs(x = 'overall score, EBCC labels (0.6 bio + 0.4 batch)', y = 'overall score, alternative labels',
       shape = 'label-informed', title = 'colon: method scores under EBCC vs best-annotator labels',
       subtitle = sprintf('best annotator by mean marker AUROC: human %s, mouse %s (a deliberately tough comparison; EBCC stays the primary labels)',
                          choice$species$human$use, choice$species$mouse$use)) +
  theme_bw() + theme(legend.position = 'bottom')
dir.create('../figures/colon/scores', showWarnings = FALSE, recursive = TRUE)
ggsave('../figures/colon/scores/label_sensitivity.png', p, width = 13, height = 9, dpi = 150)
