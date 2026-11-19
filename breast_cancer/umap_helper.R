# Figure-only UMAP of an integrated reduction, used by the 04 R scripts instead of RunUMAP.
# With only 49 shared genes, genes6 gives a kNN graph on which uwot's default spectral initialisation stalls
# or crashes (seen in 03, 04.01, 04.03, 04.04, 04.06), so the datasets in UMAP_PCA_INIT_SETS use init = 'pca'
# with otherwise RunUMAP's settings (30 neighbours, cosine, min_dist 0.3, seed 42). Every other dataset uses
# Seurat's RunUMAP unchanged. The integrations themselves are not affected: this UMAP is only drawn in figures.
UMAP_PCA_INIT_SETS = c('genes6')

run_umap = function(obj, dataset, reduction, dims, reduction.name, ...) {
  if (!dataset %in% UMAP_PCA_INIT_SETS)
    return(RunUMAP(obj, reduction = reduction, dims = dims, reduction.name = reduction.name, ...))
  set.seed(42)
  emb = uwot::umap(Embeddings(obj, reduction)[, dims, drop = FALSE], n_neighbors = 30, metric = 'cosine',
                   min_dist = 0.3, init = 'pca', n_threads = as.integer(Sys.getenv('N_THREADS', 4)))
  key = paste0(gsub('_', '', reduction.name), '_')
  rownames(emb) = colnames(obj); colnames(emb) = paste0(key, seq_len(ncol(emb)))
  obj[[reduction.name]] = CreateDimReducObject(embeddings = emb, key = key, assay = DefaultAssay(obj))
  obj
}
