# Listing A.y: aggregation into bio-conservation, batch-correction and overall scores, and
# improvement over the unintegrated data (R, dplyr/tidyr). Condensed from benchmark_figures.R.
library(dplyr); library(tidyr)

bio_metrics   = c('NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI', 'MAP')
batch_metrics = c('PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn', 'SAS')
# Graph-based methods (SAMap, BBKNN variants) have no embedding, so only the metrics they have enter
# their means: bio = NMI, ARI, cLISI; batch = iLISI, graph connectivity.

# met: one row per tissue x method x metric, column `value` (every metric on 0-1, higher = better)
score = function(met) met |>
  mutate(group = case_when(metric %in% batch_metrics ~ 'batch',
                           metric %in% bio_metrics   ~ 'bio')) |>
  filter(!is.na(group), !is.na(value)) |>
  group_by(tissue, method, group) |>
  summarise(score = mean(value), .groups = 'drop') |>            # unweighted mean within each group
  pivot_wider(names_from = group, values_from = score) |>
  mutate(overall = 0.6 * bio + 0.4 * batch)                       # scib weighting

improvement = function(scores, reference = 'Unintegrated') {
  ref = scores |> filter(method == reference) |>
    select(tissue, ref_bio = bio, ref_batch = batch, ref_overall = overall)
  scores |> filter(method != reference) |>
    left_join(ref, by = 'tissue') |>
    mutate(imp_bio     = bio     - ref_bio,                       # per tissue, minus the unintegrated
           imp_batch   = batch   - ref_batch,                     # data's score in that tissue
           imp_overall = overall - ref_overall)
}
