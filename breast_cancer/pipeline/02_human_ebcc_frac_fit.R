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

# ---- section 3b: one full EBCC fit (R = 10, all cells) for one frac, M fixed from the grid ----
# usage: Rscript pipeline/02_human_ebcc_frac_fit.R <M> <frac>
args = as.numeric(commandArgs(trailingOnly = TRUE))
M = args[1]; frac = args[2]
cat(sprintf('human: full fit with M = %s, frac = %s\n', M, frac))
fit = run_ebcc(human@meta.data, ann_cols, M = M, R = 10, beta = make_trusted_beta(M, frac))
round(fit$elbo_restarts - max(fit$elbo_restarts), 1)
out = 'data/02_annotated_data/human_ebcc_frac_fits'
dir.create(out, showWarnings = FALSE)
saveRDS(list(M = M, frac = frac, label = fit$label, posterior = fit$posterior, par = fit$par,
             elbo_restarts = fit$elbo_restarts, cells = rownames(human@meta.data)),
        file.path(out, sprintf('M%s_frac%s.rds', M, frac)))
cat('RUN FINISHED\n')
