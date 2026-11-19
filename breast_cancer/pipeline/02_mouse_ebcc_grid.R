# generated from the 02 notebook: EBCC (M, frac) grid, then the final fit with the chosen values
source('../library.R')
source('../source.R')
print('EBCC')
mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
markers_pymt = list(
  Epithelial  = c("Cdh1", "Cldn3", "Cldn4", "Cldn7", "Elf3", "Krt7", "Muc1", "Cd24a",
                  "Ehf", "Spint1", "Esrp1", "Foxa1", "Gata3", "Esr1", "Erbb2", "Erbb3"),
  Tumour      = c("Mki67", "Top2a", "Wap", "Lalba", "Csn3", "Ccnb1", "Birc5", "Pcna",
                  "Mcm2", "Cdk1"),
  Basal       = c("Krt5", "Krt14", "Krt17", "Trp63", "Acta2", "Myh11", "Cnn1", "Cdh3",
                  "Itga6", "Itgb4"),
  CAF         = c("Pdgfra", "Col1a1", "Col1a2", "Col5a1", "Col6a3", "Col14a1", "Postn",
                  "Fbln1", "Fbln2", "Dpt", "Pdpn", "Cxcl12", "Sfrp1", "Sfrp2", "Sfrp4",
                  "Mmp2", "Fn1", "Ly6a", "Pi16", "Prrx1", "Twist1", "Twist2"),
  PVL         = c("Rgs5", "Pdgfrb", "Notch3", "Cspg4", "Des", "Mcam", "Higd1b",
                  "Ndufa4l2", "Vtn"),
  Endothelial = c("Pecam1", "Cdh5", "Flt1", "Tek", "Plvap", "Vwf", "Esam", "Ptprb",
                  "Cd34", "Egfl7", "Ramp2", "Sox17", "Clec14a", "Prox1", "Lyve1", "Flt4"),
  Myeloid     = c("Adgre1", "Itgam", "Cd14", "Fcgr1", "Fcgr3", "Mrc1", "Cd163", "Itgax",
                  "C1qb", "C1qc", "Tyrobp", "Fcer1g", "Aif1", "Cx3cr1", "Ccr2", "Ly6c2",
                  "Arg1", "Spp1", "Trem2", "Apoe", "S100a8", "S100a9", "Csf3r", "Ly6g",
                  "H2-Ab1", "Cd74", "Marco", "Msr1", "Mertk"),
  T_cells     = c("Ptprc", "Cd8a", "Cd8b1", "Cd4", "Cd3g", "Cd247", "Trac", "Trbc2",
                  "Cd5", "Cd6", "Cd7", "Cd28", "Lck", "Zap70", "Itk", "Il7r", "Foxp3",
                  "Ctla4", "Pdcd1", "Gzma", "Gzmb", "Nkg7", "Tcf7", "Cd27"),
  B_cells     = c("Cd79b", "Cd22", "Ebf1", "Cr2", "Fcer2a", "Bank1", "Blk", "Cd37",
                  "Cd72", "Ighm", "Ighd", "Spib", "Tnfrsf13c"),
  Plasma      = c("Prdm1", "Tnfrsf17", "Igkc", "Ighg1", "Ighg2b", "Ighg2c", "Igha",
                  "Derl3", "Sec11c", "Fkbp11", "Txndc5", "Slamf7", "Irf4", "Pou2af1")
)

setdiff(unlist(markers_pymt), rownames(mouse))   # genes not on the panel

markers_pymt = lapply(markers_pymt, intersect, rownames(mouse))
markers_pymt = markers_pymt[lengths(markers_pymt) > 0]

mouse$celltype_seurat = mouse$celltype_seurat |>
    recode('Fibroblast' = 'CAFs')

mouse$celltype_singleR = mouse$celltype_singleR |>
    recode('Fibroblast' = 'CAFs')

ann_cols = c("celltype_seurat", "celltype_singleR", 'celltype_seurat_cross', 'celltype_singleR_cross')
# every label any annotator gives (run_ebcc needs the prior to cover all of them)
classes_all = lapply(ann_cols, function(a) as.character(mouse@meta.data[[a]])) |>
  unlist() |> na.omit() |> unique()
# hide one random annotator's label for 10% of cells
set.seed(1)
meta_ho = mouse@meta.data
ho = data.frame(cell = sample(nrow(meta_ho), round(0.1 * nrow(meta_ho))))
ho$ann   = sample(ann_cols, nrow(ho), replace = TRUE)
ho$truth = mapply(function(i, a) as.character(meta_ho[[a]][i]), ho$cell, ho$ann)
for (a in ann_cols) {
  meta_ho[[a]] = as.character(meta_ho[[a]])
  meta_ho[[a]][ho$cell[ho$ann == a]] = NA
}

# (annotator dot plots already saved)
               
trusted = list(
celltype_seurat = c('Endothelial', 'Myeloid'),
celltype_singleR  = c('Endothelial', 'CAFs', 'Myeloid'),
celltype_seurat_cross = c('PVL', 'Normal Epithelial', 'Myeloid', 'CAFs', 'Cancer Epithelial'),
celltype_singleR_cross = c('PVL', 'T-cells', 'Endothelial', 'Cancer Epithelial', 'CAFs')
)
n_k = list(
  celltype_seurat = table(mouse$celltype_seurat),
  celltype_singleR = table(mouse$celltype_singleR),
  celltype_seurat_cross = table(mouse$celltype_seurat_cross),
  celltype_singleR_cross = table(mouse$celltype_singleR_cross)
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
grid_changed = save_heldout_grid(grid, 'mouse', previous_M = 5, previous_frac = 0.5, ho = ho)
grid$grid
c(M = grid$M, frac = grid$frac, changed = grid_changed)

# ---- fit ----
# chosen by the held-out grid (see data/02_annotated_data/mouse_ebcc_heldout_grid_chosen.csv)
M = grid$M
frac = grid$frac
cat(sprintf('final fit with M = %s, frac = %s\n', M, frac))
fit = run_ebcc(mouse@meta.data, ann_cols, M = M, R = 10, beta = make_trusted_beta(M, frac))

# ---- diagnostics ----
# 1. restarts: distance of each restart's ELBO from the best (large spread = many local optima)
round(fit$elbo_restarts - max(fit$elbo_restarts), 1)

# ---- save ----
mouse$ebcc_label   = fit$label
mouse$ebcc_maxpost = rowMaxs(fit$posterior)
mouse$ebcc_entropy = -rowSums(fit$posterior * log(fit$posterior + 1e-12))
mouse$celltype     = mouse$ebcc_label
save_dotplot(mouse, markers_pymt, 'ebcc_label', 'mouse_dotplot_ebcc.png', dir = '../figures/breast_cancer/annotation')

# gene symbols + Ensembl IDs, used by 03_concatenate to map orthologs
homology_tbl = read.csv('../homology_tbl.csv')
mouse_gene_dict = homology_tbl |>
  filter(!is.na(mmusculus_homolog_orthology_type) &
           !mmusculus_homolog_orthology_type == '') |>
  dplyr::select(mmusculus_homolog_ensembl_gene,
                mmusculus_homolog_associated_gene_name) |>
  filter(!is.na(mmusculus_homolog_associated_gene_name) &
           !mmusculus_homolog_associated_gene_name == '' &
           !is.na(mmusculus_homolog_ensembl_gene) &
           !mmusculus_homolog_ensembl_gene == '') |>
  unique()
mouse[['MERFISH']]@meta.data = data.frame(gene_name = rownames(mouse)) |>
  left_join(mouse_gene_dict,
            join_by(gene_name == mmusculus_homolog_associated_gene_name))

saveRDS(fit$par,      'data/02_annotated_data/mouse_ebcc_params.rds')
saveRDS(fit$posterior, 'data/02_annotated_data/mouse_ebcc_posterior.rds')
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
as.anndata(mouse, file_path = 'data/02_annotated_data/',
           file_name = '02_mouse_annotated.h5ad', assay = 'MERFISH',
           main_layer = 'counts', other_layers = c('data'))


cat('RUN FINISHED\n')
