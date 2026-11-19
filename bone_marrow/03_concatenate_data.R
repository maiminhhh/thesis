source('../library.R')
source('../source.R')

# stages: 'inputs' (SAMap per-species h5ads), 'full', 'subset', or a removal
# fraction '0.2'/'0.4'/'0.6'/'0.8' (needs the subset stage's output). Stages can
# run as separate processes (see run_03_04.sh); no argument runs all in order.
args = commandArgs(trailingOnly = TRUE)
stages = if (length(args) > 0) args else c('inputs', 'full', 'subset', '0.2', '0.4', '0.6', '0.8')
gene_stages = setdiff(stages, c('inputs', 'full', 'subset'))

read_xenium = function(path) {
  obj = readRDS(path)
  rownames(obj) = obj[['Xenium']]@meta.data$gene_id
  obj
}
human = read_xenium('data/02_annotated_data/02_human_annotated.rds')
if (any(c('inputs', 'full') %in% stages)) {
  mouse = read_xenium('data/02_annotated_data/02_mouse_annotated.rds')
}
if (any(c('inputs', 'subset') %in% stages) || length(gene_stages) > 0) {
  mouse_subset = read_xenium('data/02_annotated_data/02_mouse_subset_annotated.rds')
}
homology_tbl = read.csv('../homology_tbl.csv')

# concatenate on one-to-one orthologues, drop empty cells, embed and (if
# resolutions are given) cluster
make_concat = function(human, mouse, resolution = NULL) {
  concat = concatenate(human, mouse, homology_tbl)
  VariableFeatures(concat) = rownames(concat)
  concat = JoinLayers(concat)
  concat$nCount_Xenium = colSums(x = concat[['Xenium']]$counts)
  concat = concat |>
    subset(nCount_Xenium > 0) |>
    RunPCA() |>
    RunUMAP(dims = 1:20, uwot.init = Sys.getenv('UMAP_INIT', 'spectral'))   # UMAP_INIT=random if spectral hangs (genes2)
  if (!is.null(resolution)) {
    concat = concat |>
      FindNeighbors(dims = 1:20) |>
      FindClusters(resolution = resolution, algorithm = 4, random.seed = 1)
  }
  concat
}

save_both = function(obj, dir, name) {
  SaveSeuratRds(obj, file.path(dir, paste0(name, '.rds')))
  as.anndata(obj, file_path = dir, file_name = paste0(name, '.h5ad'),
             assay = 'Xenium', main_layer = 'counts',
             other_layers = c('data', 'scale.data'))
}

# per-species inputs for SAMap on the full/dissected data, keyed by Ensembl ID
# like the genes* folders (the 02_*.h5ad files use gene symbols)
if ('inputs' %in% stages) {
  for (x in list(list(human, '03_human.h5ad'), list(mouse, '03_mouse_full.h5ad'),
                 list(mouse_subset, '03_mouse_subset.h5ad'))) {
    as.anndata(x[[1]], file_path = 'data/03_concatenated_data/', file_name = x[[2]],
               assay = 'Xenium', main_layer = 'counts', other_layers = c('data', 'scale.data'))
  }
}

# full section
# ElbowPlot(concat, ndims = 50)
# ggsave('elbow_plot.png', ElbowPlot(concat, ndims = 50))
if ('full' %in% stages) {
  concat = make_concat(human, mouse, resolution = c(0.4, 0.6, 0.8, 1.2))
  save_both(concat, 'data/03_concatenated_data/', '03_concatenated_data_full')
}

# subsection
if ('subset' %in% stages) {
  concat = make_concat(human, mouse_subset, resolution = c(0.2, 0.4, 0.6, 0.8, 1.2))
  save_both(concat, 'data/03_concatenated_data/', '03_concatenated_data_subset')
}

# subset genes
# remove 20/40/60/80% of the shared genes (saved to genes2/4/6/8). Removal is
# nested: one random order of the shared genes, each level drops a longer prefix.
# Clusters are not computed here: nothing downstream uses them.
if (length(gene_stages) > 0) {
  if (!'subset' %in% stages) {
    concat = readRDS('data/03_concatenated_data/03_concatenated_data_subset.rds')
  }
  shared_genes = concat[['Xenium']]@meta.data
  stopifnot(identical(shared_genes$hsapiens_homolog_ensembl_gene, rownames(concat)))
  set.seed(1)
  removal_order = sample(nrow(shared_genes))
}

# counts are unchanged; data/scale.data must be recomputed on the genes that remain
renormalize = function(obj) {
  obj$nCount_Xenium = colSums(x = obj[['Xenium']]$counts)
  obj = obj |>
    subset(nCount_Xenium > 0)
  VariableFeatures(obj) = rownames(obj)
  obj |>
    NormalizeData() |>
    ScaleData()
}

for (p in as.numeric(gene_stages)) {
  out_dir = paste0('data/genes', p * 10, '/')
  removed = shared_genes[removal_order[seq_len(round(p * nrow(shared_genes)))], ]
  human_p = human[!rownames(human) %in% removed$hsapiens_homolog_ensembl_gene, ] |>
    renormalize()
  mouse_p = mouse_subset[!rownames(mouse_subset) %in% removed$mmusculus_homolog_ensembl_gene, ] |>
    renormalize()
  concat_p = make_concat(human_p, mouse_p)
  print(str_glue('{out_dir}: removed {nrow(removed)} of {nrow(shared_genes)} shared genes, ',
                 '{nrow(concat_p)} remain'))

  save_both(human_p, out_dir, '02_human_annotated')
  save_both(mouse_p, out_dir, '02_mouse_annotated')
  save_both(concat_p, out_dir, '03_concatenated_data')
}
