# generated: section 3b of ../EBCC_PIPELINE.md for mouse. M from the held-out grid (1-SE rule on the frac = 0
# column), frac by marker AUROC on a 50/50 split; the saved full fit for the chosen frac becomes the final fit.
source('../library.R')
source('../source.R')
source('ebcc_grid.R')
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


# EBCC label -> marker group of `markers_pymt` (the marker lists themselves are unchanged)
label_map = c('B-cells' = 'B_cells', 'CAFs' = 'CAF', 'Cancer Epithelial' = 'Tumour', 'Endothelial' = 'Endothelial',
  'Myeloid' = 'Myeloid', 'Normal Epithelial' = 'Basal', 'Plasmablasts' = 'Plasma', 'PVL' = 'PVL', 'T-cells' = 'T_cells')

grid = read.csv('data/02_annotated_data/mouse_ebcc_heldout_grid.csv')
M = choose_M_frac0(grid)
fracs = c(0, 0.25, 0.5, 1, 2, 5)
fits = lapply(fracs, function(f) readRDS(sprintf('data/02_annotated_data/mouse_ebcc_frac_fits/M%s_frac%s.rds', M, f)))
names(fits) = as.character(fracs)
stopifnot(all(vapply(fits, function(f) identical(f$cells, colnames(mouse)), logical(1))))

X = LayerData(mouse, assay = DefaultAssay(mouse), layer = 'data')
res = frac_auroc(X, fits, markers_pymt, label_map)
res$summary
changed = save_frac_auroc(res, 'mouse', M = M, previous_frac = 0.5, previous_M = 5)
cat(sprintf('mouse: M = %s, frac = %s (selection AUROC %.4f, held-out %.4f), changed = %s\n', M, res$frac,
            res$summary$mean_auroc_selection[res$summary$chosen], res$summary$mean_auroc_heldout[res$summary$chosen], changed))

# final fit = the saved fit for the chosen frac (no refit)
frac = res$frac
fit = fits[[as.character(frac)]]
cat(sprintf('final fit: M = %s, frac = %s\n', M, frac))
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
