# generated: section 3b of ../EBCC_PIPELINE.md for human. M from the held-out grid (1-SE rule on the frac = 0
# column), frac by marker AUROC on a 50/50 split; the saved full fit for the chosen frac becomes the final fit.
source('../library.R')
source('../source.R')
source('ebcc_grid.R')
human = readRDS('data/02_annotated_data/02_human_annotated.rds')

markers = list(
  Epithelial   = c("EPCAM", "CDH1", "KRT7", "ELF3", "CLDN4"),
  Basal        = c("KRT5", "KRT14", "TP63", "ACTG2"),
  Tumour       = c("ERBB2", "MKI67", "TOP2A", "CEACAM6", "S100P"),
  CAF          = c("PDGFRA", "FAP", "POSTN", "FN1", "FBLN1"),
  PVL          = c("RGS5", "PDGFRB", "MCAM", "MYH11"),
  Endothelial  = c("PECAM1", "VWF", "CDH5", "CLDN5", "KDR"),
  T_NK         = c("PTPRC", "CD3E", "CD2", "NKG7", "GNLY"),
  B            = c("MS4A1", "CD19", "CD79A", "PAX5"),
  Plasmablast  = c("MZB1", "DERL3", "TNFRSF17", "IGHG3"),
  Myeloid      = c("CD68", "CD163", "CD14", "LYZ", "C1QA")
)

setdiff(unlist(markers), rownames(human))   # genes not on the panel

markers = lapply(markers, intersect, rownames(human))
markers = markers[lengths(markers) > 0]

# EBCC label -> marker group of `markers` (the marker lists themselves are unchanged)
label_map = c('B-cells' = 'B', 'CAFs' = 'CAF', 'Cancer Epithelial' = 'Tumour', 'Endothelial' = 'Endothelial',
  'Myeloid' = 'Myeloid', 'Normal Epithelial' = 'Basal', 'Plasmablasts' = 'Plasmablast', 'PVL' = 'PVL', 'T-cells' = 'T_NK')

grid = read.csv('data/02_annotated_data/human_ebcc_heldout_grid.csv')
M = choose_M_frac0(grid)
fracs = c(0, 0.25, 0.5, 1, 2, 5)
fits = lapply(fracs, function(f) readRDS(sprintf('data/02_annotated_data/human_ebcc_frac_fits/M%s_frac%s.rds', M, f)))
names(fits) = as.character(fracs)
stopifnot(all(vapply(fits, function(f) identical(f$cells, colnames(human)), logical(1))))

X = LayerData(human, assay = DefaultAssay(human), layer = 'data')
res = frac_auroc(X, fits, markers, label_map)
res$summary
changed = save_frac_auroc(res, 'human', M = M, previous_frac = 0, previous_M = 3)
cat(sprintf('human: M = %s, frac = %s (selection AUROC %.4f, held-out %.4f), changed = %s\n', M, res$frac,
            res$summary$mean_auroc_selection[res$summary$chosen], res$summary$mean_auroc_heldout[res$summary$chosen], changed))
if (!changed) { cat('human: M and frac unchanged, the saved final fit stands\nRUN FINISHED\n'); quit(save = 'no') }

# final fit = the saved fit for the chosen frac (no refit)
frac = res$frac
fit = fits[[as.character(frac)]]
cat(sprintf('final fit: M = %s, frac = %s\n', M, frac))
# ---- diagnostics ----
# restarts: distance of each restart's ELBO from the best (large spread = many local optima)
round(fit$elbo_restarts - max(fit$elbo_restarts), 1)

# ---- save ----
human$ebcc_label   = fit$label
human$ebcc_maxpost = rowMaxs(fit$posterior)
human$ebcc_entropy = -rowSums(fit$posterior * log(fit$posterior + 1e-12))
save_dotplot(human, markers, 'ebcc_label', 'human_dotplot_ebcc.png', dir = '../figures/breast_cancer/annotation')

saveRDS(fit$par,       'data/02_annotated_data/human_ebcc_params.rds')
saveRDS(fit$posterior, 'data/02_annotated_data/human_ebcc_posterior.rds')
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')


## -----------------------------------------------------------------------------
human$celltype_seurat |> table(human$ebcc_label)
human$celltype_seurat_cross |> table(human$ebcc_label)
ImageDimPlot(human, group.by = 'ebcc_label')
DimPlot(human, reduction = 'umap', group.by = 'ebcc_label', shuffle = TRUE)


## -----------------------------------------------------------------------------
human$celltype = human$ebcc_label


## -----------------------------------------------------------------------------
homology_tbl = read.csv('../homology_tbl.csv')
human_gene_dict = homology_tbl |> 
  arrange(mmusculus_homolog_ensembl_gene) |> 
  filter(!is.na(mmusculus_homolog_orthology_type) & 
           !mmusculus_homolog_orthology_type == '') |> 
  dplyr::select(ensembl_gene_id, external_gene_name) |> 
  filter(!is.na(external_gene_name)) |> 
  filter(!external_gene_name == '') |> 
  unique()

# rebuilt from the gene names, so rerunning this chunk never duplicates columns
human[['MERFISH']]@meta.data = data.frame(gene_name = rownames(human)) |> 
  left_join(human_gene_dict, join_by(gene_name == external_gene_name))
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')


## -----------------------------------------------------------------------------
write_h5ad(human, 'data/02_annotated_data/02_human_annotated.h5ad', mode = 'w')


## -----------------------------------------------------------------------------
ImageDimPlot(human, group.by = 'celltype', axes = TRUE)


cat('RUN FINISHED\n')
