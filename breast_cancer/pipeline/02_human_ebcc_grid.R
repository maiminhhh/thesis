# generated from the 02 notebook: EBCC (M, frac) grid, then the final fit with the chosen values
source('../library.R')
source('../source.R')
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

ann_cols = c('celltype_seurat', 'celltype_singleR', 'celltype_seurat_cross', 'celltype_singleR_cross')
# every label any annotator gives (run_ebcc needs the prior to cover all of them)
classes_all = lapply(ann_cols, function(a) as.character(human@meta.data[[a]])) |> 
  unlist() |> na.omit() |> unique()

# hide one random annotator's label for 10% of cells
set.seed(1)
meta_ho = human@meta.data
ho = data.frame(cell = sample(nrow(meta_ho), round(0.1 * nrow(meta_ho))))
ho$ann   = sample(ann_cols, nrow(ho), replace = TRUE)
ho$truth = mapply(function(i, a) as.character(meta_ho[[a]][i]), ho$cell, ho$ann)
for (a in ann_cols) {
  meta_ho[[a]] = as.character(meta_ho[[a]])
  meta_ho[[a]][ho$cell[ho$ann == a]] = NA
}

trusted = list(
celltype_seurat = c('T-cells', 'PVL' ,'Plasmablasts', 'Normal Epithelial', 'Myeloid', 'Cancer Epithelial', 'CAFs', 'B-cells'),
celltype_singleR  = c('T-cells', 'Normal Epithelial', 'PVL', 'Plasmablasts', 'Myeloid', 'Cancer Epithelial', 'CAFs'),
celltype_seurat_cross = c('Myeloid'),
celltype_singleR_cross = c('T-cells', 'Myeloid')
)
n_k = list(
  celltype_seurat = table(human$celltype_seurat),
  celltype_singleR = table(human$celltype_singleR),
  celltype_seurat_cross = table(human$celltype_seurat_cross),
  celltype_singleR_cross = table(human$celltype_singleR_cross)
)

# default priors, strengthened on the diagonal for each annotator's trusted classes
# (trusted classes the annotator never called are skipped)
make_trusted_beta = function(M, frac) {
  beta = make_beta_prior(classes_all, ann_cols)
  for (a in names(trusted)) for (k in trusted[[a]]) {
    n = if (k %in% names(n_k[[a]])) n_k[[a]][[k]] else 0
    if (n > 0) beta[[a]][k, k] = 4 + frac * n / M
  }
  beta
}

# fit with M subtypes, then score how well the hidden labels are predicted
heldout_ll = function(M, frac, R = 3) {
  f = run_ebcc(meta_ho, ann_cols, M = M, R = R, keep_rho = TRUE, beta = make_trusted_beta(M, frac))
  lp = unlist(lapply(ann_cols, function(a) {
    s = ho[ho$ann == a, ]
    l = match(s$truth, f$classes); ok = !is.na(l)
    # P(hidden label) = sum over class-subtypes of P(cell in it) x P(annotator gives that label)
    log(rowSums(f$rho[s$cell[ok], , drop = FALSE] *
                t(f$par$v[[a]][, l[ok], drop = FALSE])))
  }))
  lp
}

# ---- joint (M, frac) held-out grid ----
# M in 1, 2, 3, 5, 10 x frac in 0, 0.25, 0.5, 1, 2; 3 restarts per fit; one process per frac.
# Rule: best mean held-out log-likelihood; among (M, frac) within 1 SE of it, smallest M then smallest frac.
source('ebcc_grid.R')
grid = heldout_grid(heldout_ll, Ms = c(1, 2, 3, 5, 10), fracs = c(0, 0.25, 0.5, 1, 2), threads = 4)
grid_changed = save_heldout_grid(grid, 'human', previous_M = 5, previous_frac = 0.5, ho = ho)
grid$grid
c(M = grid$M, frac = grid$frac, changed = grid_changed)
if (!grid_changed) { cat('human: grid kept M = 5, frac = 0.5, so the existing final fit stands\nRUN FINISHED\n'); quit(save = 'no') }

# ---- fit ----
# chosen by the held-out grid (see data/02_annotated_data/human_ebcc_heldout_grid_chosen.csv)
M = grid$M
frac = grid$frac
cat(sprintf('final fit with M = %s, frac = %s\n', M, frac))
fit = run_ebcc(human@meta.data, ann_cols, M = M, R = 10, beta = make_trusted_beta(M, frac))

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
