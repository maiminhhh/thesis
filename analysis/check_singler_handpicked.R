# How many genes SingleR uses to identify each cell type, and how many of the hand-picked markers are among them.
# SingleR (de.method = 't', default de.n = 10) picks, for every ordered pair of types (A, B), the 10 genes most
# up-regulated in A vs B. A type's genes = union over its comparisons (analysis/check_singler_genes.R, figures/cross_tissue/
# singler_genes/). All types are first scored by Spearman correlation over the union of every type's genes (the
# shared set); fine-tuning then re-correlates using the pairwise genes among the closest types.
# Hand-picked markers: annotation_markers (the lists used for marker AUROC), matched to SingleR types ignoring case and
# punctuation. Same-species reference; breast cancer mouse left out. Reads saved CSVs and marker scripts only.
# Usage: Rscript45 analysis/check_singler_handpicked.R -> figures/cross_tissue/singler_handpicked.csv, printed summary
ex = parse('analysis/vote_vs_ebcc.R')
stop_at = which(startsWith(sapply(ex, \(e) deparse(e)[1]), 'args = commandArgs'))
for (e in ex[seq_len(stop_at - 1)]) eval(e, globalenv())

sg = list.files('figures/cross_tissue/singler_genes', '\\.csv$', full.names = TRUE) |>
  lapply(read_csv, show_col_types = FALSE) |> bind_rows()
panel = read_csv('figures/cross_tissue/raw_gene_detection.csv', show_col_types = FALSE) |> filter(feature_type == 'Gene Expression')
res = sg |> group_by(tissue, species) |> group_modify(function(g, k) {
  sr = setNames(strsplit(g$genes, ';'), g$type)                    # SingleR's genes per type (symbols)
  shared = unique(unlist(sr))
  hp = annotation_markers(tissues[[k$tissue]]$annotation[[k$species]])
  hp = setNames(hp, key(names(hp)))
  pg = panel$feature[panel$tissue == k$tissue & panel$species == k$species]
  tibble(type = names(sr), n_singler = lengths(sr), n_shared_set = length(shared),
         n_handpicked = vapply(key(names(sr)), \(t) length(unique(hp[[t]])), integer(1)),
         n_handpicked_in_panel = vapply(key(names(sr)), \(t) length(intersect(unique(hp[[t]]), pg)), integer(1)),
         n_in_type_set = vapply(seq_along(sr), \(i) length(intersect(unique(hp[[key(names(sr))[i]]]), sr[[i]])), integer(1)),
         n_in_shared_set = vapply(key(names(sr)), \(t) length(intersect(unique(hp[[t]]), shared)), integer(1)),
         in_type_set = vapply(seq_along(sr), \(i) paste(intersect(unique(hp[[key(names(sr))[i]]]), sr[[i]]), collapse = ';'), character(1)),
         has_handpicked = key(names(sr)) %in% names(hp))
}) |> ungroup() |> mutate(tissue = factor(tissue, levels = tissue_order)) |> arrange(tissue, species, type)
write_csv(res, 'figures/cross_tissue/singler_handpicked.csv')
print(as.data.frame(res |> select(tissue, species, type, n_singler, n_handpicked, n_handpicked_in_panel, n_in_type_set, n_in_shared_set)), right = FALSE)
print(as.data.frame(res |> filter(has_handpicked) |> group_by(tissue, species) |>
  summarise(shared_set = first(n_shared_set), handpicked_in_panel = sum(n_handpicked_in_panel),
            in_type_set = sum(n_in_type_set), in_shared_set = sum(n_in_shared_set), .groups = 'drop')))
cat('types without hand-picked markers:', paste(res$tissue[!res$has_handpicked], res$species[!res$has_handpicked], res$type[!res$has_handpicked]), '\n')
