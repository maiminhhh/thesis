# Does sparsity drive the concentration of labels into one dominant type? Within each Seurat label (so cell type is
# held roughly fixed), cells are split into quintiles of panel genes detected (count > 0); per quintile, the share
# that SingleR puts into its dataset-wide dominant label, and the same for Seurat's dominant label among SingleR
# labels (the reverse direction). Prediction if sparsity drives SingleR's concentration: the SingleR share rises as
# genes detected fall, more than Seurat's. Same-species reference; SingleR pruned calls left out.
# Run from the repository root. Reads saved objects only.
# Usage: Rscript analysis/check_sparsity_concentration.R <tissue> <species> -> figures/cross_tissue/sparsity_concentration/<tissue>_<species>.csv
#        Rscript analysis/check_sparsity_concentration.R                    -> figures/cross_tissue/sparsity_concentration.csv, printed
ex = parse('analysis/vote_vs_ebcc.R')
stop_at = which(startsWith(sapply(ex, \(e) deparse(e)[1]), 'args = commandArgs'))
for (e in ex[seq_len(stop_at - 1)]) eval(e, globalenv())

out_dir = 'figures/cross_tissue/sparsity_concentration'
args = commandArgs(TRUE)
if (length(args) == 2) {
  dir.create(out_dir, showWarnings = FALSE)
  a = tissues[[args[1]]]$annotation[[args[2]]]
  obj = readRDS(a$obj)
  X = LayerData(obj, assay = DefaultAssay(obj), layer = 'counts')
  m = obj@meta.data
  rm(obj); gc()
  d = tibble(genes = Matrix::colSums(X > 0), seurat = as.character(m$celltype_seurat),
             singleR = as.character(m$celltype_singleR), delta = m$labeldelta_singleR, score = m$labelscore_seurat) |>
    filter(!is.na(seurat), !is.na(singleR))
  dom_singleR = names(which.max(table(d$singleR))); dom_seurat = names(which.max(table(d$seurat)))
  # quintile of genes detected within each label (ties broken at random so quintiles have equal size)
  set.seed(1)
  q_within = \(g, lab) ave(g + runif(length(g), 0, 1e-3), lab, FUN = \(x) dplyr::ntile(x, 5))
  one = function(lab_col, other_col, dom, direction) {
    d |> filter(.data[[lab_col]] != dom) |>                              # cells the first annotator does NOT call dominant
      group_by(lab = .data[[lab_col]]) |> filter(n() >= 200) |> ungroup() |>
      mutate(q = q_within(genes, lab)) |>
      group_by(q) |> summarise(n = n(), median_genes = median(genes), to_dominant = mean(.data[[other_col]] == dom),
                               agree = mean(.data[[lab_col]] == .data[[other_col]]), .groups = 'drop') |>
      mutate(direction = direction, dominant = dom)
  }
  res = bind_rows(one('seurat', 'singleR', dom_singleR, 'SingleR pulls into its dominant label'),
                  one('singleR', 'seurat', dom_seurat, 'Seurat pulls into its dominant label')) |>
    mutate(tissue = args[1], species = args[2], .before = 1)
  write_csv(res, file.path(out_dir, sprintf('%s_%s.csv', args[1], args[2])))
  quit(save = 'no')
}
d = purrr::map_dfr(list.files(out_dir, '\\.csv$', full.names = TRUE), read_csv, show_col_types = FALSE)
write_csv(d, 'figures/cross_tissue/sparsity_concentration.csv')
print(as.data.frame(d |> select(tissue, species, direction, dominant, q, median_genes, to_dominant, agree)), digits = 2)
