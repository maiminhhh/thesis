# Best single annotator per dataset and species, by mean marker AUROC on the SELECTION half of cells
# (the half used to choose f; the complement of the evaluation half in analysis/vote_vs_ebcc.R). Candidates: the four
# annotators only (no vote, no EBCC); each scored on the cells it labels, mean over the types it calls.
# Also records which confidence score each annotator has saved (Seurat: labelscore_*, SingleR: labeldelta_*).
# Usage:  Rscript analysis/best_annotator_select.R <tissue> <species>  -> figures/cross_tissue/best_annotator_select/<tissue>_<species>.csv
#         Rscript analysis/best_annotator_select.R                     -> figures/cross_tissue/best_annotator_select.csv

ex = parse('analysis/vote_vs_ebcc.R')
stop_at = which(startsWith(sapply(ex, \(e) deparse(e)[1]), 'args = commandArgs'))
for (e in ex[seq_len(stop_at - 1)]) eval(e, globalenv())

out_dir = 'figures/cross_tissue/best_annotator_select'
ann_names = c(seurat = 'Seurat', singleR = 'SingleR', seurat_cross = 'Seurat (cross)', singleR_cross = 'SingleR (cross)')
score_cols = c(seurat = 'labelscore_seurat', singleR = 'labeldelta_singleR',
               seurat_cross = 'labelscore_seurat_cross', singleR_cross = 'labeldelta_singleR_cross')

args = commandArgs(TRUE)
if (length(args) == 2) {
  tissue = args[1]; species = args[2]
  a = tissues[[tissue]]$annotation[[species]]
  obj = readRDS(a$obj)
  assay = DefaultAssay(obj)
  rownames(obj) = obj[[assay]]@meta.data$gene_name
  X = LayerData(obj, assay = assay, layer = 'data')
  meta = obj@meta.data
  rm(obj); gc()
  set.seed(1)
  select = seq_len(ncol(X)) %in% sample(ncol(X), floor(ncol(X) / 2))   # same split as analysis/vote_vs_ebcc.R
  markers = lapply(annotation_markers(a), intersect, rownames(X))
  markers = markers[lengths(markers) > 0]
  score = sapply(markers, \(g) Matrix::colMeans(X[g, select, drop = FALSE]))
  colnames(score) = key(colnames(score))
  rm(X); gc()
  res = purrr::map_dfr(names(ann_names), function(n) {
    l = key(as.character(meta[[annotators[[n]]]])[select])
    idx = which(!is.na(l))
    types = colnames(score)[vapply(colnames(score), \(k) { p = l[idx] == k; any(p) && !all(p) }, logical(1))]
    sc = score_cols[[n]]
    tibble(tissue = tissue, species = species, annotator = ann_names[[n]], n_select = sum(select),
           n_labelled = length(idx), n_types = length(types),
           mean_auroc = mean(vapply(types, \(k) auroc(score[idx, k], l[idx] == k), numeric(1))),
           score_col = sc, score_saved = sc %in% colnames(meta),
           score_na = if (sc %in% colnames(meta)) mean(is.na(meta[[sc]])) else NA_real_)
  })
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(res, file.path(out_dir, sprintf('%s_%s.csv', tissue, species)))
  quit(save = 'no')
}

res = map_dfr(tissue_order, \(t) map_dfr(c('human', 'mouse'), \(sp)
  read_csv(file.path(out_dir, sprintf('%s_%s.csv', t, sp)), show_col_types = FALSE))) |>
  group_by(tissue, species) |> mutate(best = mean_auroc == max(mean_auroc)) |> ungroup()
write_csv(res, 'figures/cross_tissue/best_annotator_select.csv')
print(as.data.frame(res |> select(tissue, species, annotator, mean_auroc, best, score_saved, score_na)))
