# generated from 03_concatenate.qmd: plots and plot-only PCA/UMAP of the species objects removed
source('../library.R')
source('../source.R')
homology_tbl = read.csv('../homology_tbl.csv')
# remove 20/40/60/80% of the shared genes from the dissected dataset (saved to genes2/4/6/8;
# genes8 = 80% removed). Removal is nested: one random order of the shared genes, each level
# drops a longer prefix. Each species is renormalised on the genes that remain, then
# concatenated again (as in bone_marrow, see ../EBCC_PIPELINE.md section 6).
# the 02 PCA/UMAP of the species objects are not used here (03 makes its own) and break gene / cell renames
drop_reductions = function(obj) { for (r in Reductions(obj)) obj[[r]] = NULL; obj }
human = readRDS('data/dissected_mouse/02_human_annotated.rds') |> drop_reductions()
mouse = readRDS('data/dissected_mouse/02_mouse_annotated.rds') |> drop_reductions()
rownames(human) = coalesce(human[['MERFISH']]@meta.data$ensembl_gene_id,
                           human[['MERFISH']]@meta.data$gene_name)
rownames(mouse) = coalesce(mouse[['MERFISH']]@meta.data$mmusculus_homolog_ensembl_gene,
                           mouse[['MERFISH']]@meta.data$gene_name)
shared_genes = readRDS('data/dissected_mouse/03_concatenated_data.rds')[['MERFISH']]@meta.data
set.seed(1)
removal_order = sample(nrow(shared_genes))

# counts are unchanged; data/scale.data are recomputed on the genes that remain
renormalize = function(obj) {
  obj$nCount_MERFISH = colSums(x = obj[['MERFISH']]$counts)
  obj = obj |>
    subset(nCount_MERFISH > 0)
  VariableFeatures(obj) = rownames(obj)
  obj |>
    NormalizeData() |>
    ScaleData()
}

# species objects carry the hsapiens_/mmusculus_ prefixes; concatenate() adds them itself
make_concat = function(human, mouse) {
  human = RenameCells(human, new.names = sub('^hsapiens_', '', colnames(human)))
  mouse = RenameCells(mouse, new.names = sub('^mmusculus_', '', colnames(mouse)))
  concat = concatenate(human, mouse, homology_tbl, tech_human = 'MERFISH', tech_mouse = 'MERFISH')
  VariableFeatures(concat) = rownames(concat)
  concat = JoinLayers(concat)
  concat$nCount_MERFISH = colSums(x = concat[['MERFISH']]$counts)
  concat = concat |>
    subset(nCount_MERFISH > 0)
  # PCA + UMAP of the unintegrated data, used by the figures
  concat |>
    RunPCA(npcs = min(50, nrow(concat) - 1)) |>
    RunUMAP(dims = 1:20)
}

save_dataset = function(x, out) {
  dir.create(out, showWarnings = FALSE)
  SaveSeuratRds(x$mouse, file.path(out, '02_mouse_annotated.rds'))
  SaveSeuratRds(x$human, file.path(out, '02_human_annotated.rds'))
  SaveSeuratRds(x$concat, file.path(out, '03_concatenated_data.rds'))
  as.anndata(x$mouse, file_path = paste0(out, '/'), file_name = '02_mouse_annotated.h5ad', assay = 'MERFISH',
             main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')
  as.anndata(x$human, file_path = paste0(out, '/'), file_name = '02_human_annotated.h5ad', assay = 'MERFISH',
             main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')
  as.anndata(x$concat, file_path = paste0(out, '/'), file_name = '03_concatenated_data.h5ad', assay = 'MERFISH',
             main_layer = 'counts', other_layers = c('data', 'scale.data'))
}

# p = fraction of the shared genes removed
reduce_genes = function(p) {
  out = file.path('data', paste0('genes', round(p * 10)))
  removed = shared_genes[removal_order[seq_len(round(p * nrow(shared_genes)))], ]
  human_p = human[!rownames(human) %in% removed$hsapiens_homolog_ensembl_gene, ] |> renormalize()
  mouse_p = mouse[!rownames(mouse) %in% removed$mmusculus_homolog_ensembl_gene, ] |> renormalize()
  concat_p = make_concat(human_p, mouse_p)
  # species objects keep the cells of the concatenated data, as for full_mouse / dissected_mouse
  human_p = human_p[, intersect(colnames(human_p), colnames(concat_p))]
  mouse_p = mouse_p[, intersect(colnames(mouse_p), colnames(concat_p))]
  print(str_glue('{out}: removed {nrow(removed)} of {nrow(shared_genes)} shared genes, ',
                 '{nrow(concat_p)} remain; {ncol(concat_p)} cells'))
  save_dataset(list(human = human_p, mouse = mouse_p, concat = concat_p), out)
}

# one gene-removal level per process: Rscript 03_partB_reduce_genes.R 0.2
reduce_genes(as.numeric(commandArgs(trailingOnly = TRUE)[1]))

