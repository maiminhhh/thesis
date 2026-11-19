source('../library.R')
source('../source.R')
mouse = readRDS('data/01_processed_data/01_mouse_processed.rds')
mouse = mouse  |> 
  NormalizeData() |> 
  ScaleData()
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
print('Done setting up!')

# within species reference
# ref source: https://altanalyze.org/MarrowAtlas/
mouse_ref = read_h5ad('data/02_annotation_refs/adata_combined_60k_rna_adt.h5ad')
mouse_ref = mouse_ref$as_Seurat()
# X is log-normalized (no raw counts in this file); LogNormalize is per-cell scale
# invariant, so expm1(X) stands in for counts. Drop the ADT features first.
mouse_ref = mouse_ref[!grepl('^ADT-', rownames(mouse_ref)),]
mouse_ref[['RNA']]$counts = expm1(mouse_ref[['RNA']]$X)
mouse_ref[['RNA']]$X = NULL
mouse_ref = mouse_ref |>
    NormalizeData()

mouse_l1 = mouse_ref$`Level 1` |> 
    as.character()
mouse_ref$celltype = case_when(
    mouse_l1 %in% c('HSPC', 'HSCP', 'MultiLin') ~ 'HSPC',
    mouse_l1 == 'MkP' ~ 'MkP',
    mouse_l1 == 'ERP' ~ 'Erythroid',
    # Mast/Mac folded into classes the human reference also has, so both references share one label set
    mouse_l1 %in% c("preNeu", "Baso", "Eosinophil", "Mast") ~ "Granulocyte",
    mouse_l1 %in% c('Mono', 'Mac') ~ 'Monocyte',
    mouse_l1 == 'DC' ~ 'DC',
    mouse_l1 == 'B cell' ~ 'B cell',
    mouse_l1 == 'T cell' ~ 'T cell',
    mouse_l1 == 'ILC' ~ 'NK/ILC'
)  |> factor()

SaveSeuratRds(mouse_ref, 'data/02_annotation_refs/02_mouse_ref.rds')
print('Done pre-processing reference!')

mouse_intersect = intersect(rownames(mouse), rownames(mouse_ref))
mouse_annot = mouse[mouse_intersect,]
mouse_ref = mouse_ref[mouse_intersect,]

# label transfer
mouse_seurat_label = LabelTransfer(data = mouse_annot, ref = mouse_ref,
                                   celltype_column = 'celltype')
mouse$celltype_seurat = mouse_seurat_label$cell_type_transfer |> as.factor()
mouse$labelscore_seurat = mouse_seurat_label$score_transfer
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
print('Done Label Transfer!')

# singleR
mouse_singleR_label = LabelSingleR(data = mouse_annot, ref = mouse_ref, 
                                   celltype_column = 'celltype',
                                   tech = 'Xenium') 
mouse$celltype_singleR = mouse_singleR_label$cell_type_singleR |>
    as.factor()
mouse$labeldelta_singleR = mouse_singleR_label$delta_singleR
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
print('Done SingleR!')

# cross-species reference (HUMAN reference, mapped onto mouse genes)
mouse_ref = readRDS('data/02_annotation_refs/02_human_ref.rds')
mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
homology_tbl = read.csv('../homology_tbl.csv')

mouse_gene_dict = homology_tbl |> 
    filter(!is.na(ensembl_gene_id))  |> 
    filter(!is.na(external_gene_name)) |>
    dplyr::select(ensembl_gene_id, external_gene_name) |>
    filter(!ensembl_gene_id == '')  |> 
    filter(!external_gene_name == '')  |> 
    distinct()

gene_name = rownames(mouse_ref)
ens = mouse_gene_dict$ensembl_gene_id[
    match(gene_name,mouse_gene_dict$external_gene_name)
]
gene_id = coalesce(
    ens,
    gene_name
)
rownames(mouse_ref) = gene_id

rownames(mouse) = mouse[['Xenium']]@meta.data$gene_id
mouse_ref_mouse = concatenate(mouse_ref, mouse, homology_tbl)
rownames(mouse_ref_mouse) = mouse_ref_mouse[['Xenium']]@meta.data$mmusculus_homolog_ensembl_gene
DefaultAssay(mouse_ref_mouse) = 'RNA'
DefaultLayer(mouse_ref_mouse[['RNA']]) = 'data'
mouse_ref_mouse = mouse_ref_mouse |> 
    subset(orig.ident == 'hsapiens')
mouse_ref_mouse$nCount_RNA = colSums(x = mouse_ref_mouse[['RNA']]$counts)
mouse_ref_mouse = mouse_ref_mouse |> 
    subset(nCount_RNA > 0)
SaveSeuratRds(mouse_ref_mouse, 'data/02_annotation_refs/human_processed_mouse_use.rds')
print('Done pre-processing cross-species reference!')

# label transfer
mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
rownames(mouse) = mouse[['Xenium']]@meta.data$gene_id
mouse_ref_mouse = readRDS('data/02_annotation_refs/human_processed_mouse_use.rds')
mouse_intersect = intersect(rownames(mouse), rownames(mouse_ref_mouse))
mouse_annot = mouse[mouse_intersect,]
xref = mouse_ref_mouse[mouse_intersect,]

mouse_seurat_label = LabelTransfer(data = mouse_annot, ref = xref,
                                   celltype_column = 'celltype')
mouse$celltype_seurat_cross = mouse_seurat_label$cell_type_transfer |> as.factor()
mouse$labelscore_seurat_cross = mouse_seurat_label$score_transfer
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
print('Done Label Transfer - Cross-species Reference!')

# singleR
mouse_singleR_label = LabelSingleR(data = mouse_annot, ref = xref,
                                   celltype_column = 'celltype',
                                   tech = 'Xenium')
mouse$celltype_singleR_cross = mouse_singleR_label$cell_type_singleR |>
    as.factor()
mouse$labeldelta_singleR_cross = mouse_singleR_label$delta_singleR
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
print('Done SingleR - Cross-species Reference!')

print('EBCC')
mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
rownames(mouse) = mouse[['Xenium']]@meta.data$gene_name

# mouse bone marrow markers chosen from the Xenium panel. Macrophage, Endothelial
# and Chondro/Osteo are not annotation classes: they show where those cells end up.
markers = list(
  HSPC            = c("Cd34", "Kit", "Ly6a", "Slamf1"),
  MkP             = c("Pf4", "Ppbp", "Vwf"),
  Erythroid       = c("Gata1", "Epor", "Klf1", "Hemgn", "Tfrc", "Gypa", "Slc4a1", "Hbb-bt", "Rhd"),
  Granulocyte     = c("Mpo", "Elane", "Ctsg", "Ly6g", "S100a9", "Mmp8", "Mmp9", "Csf3r", "Cxcr2", "Fcnb",
                      "Prg2", "Cpa3"),
  Monocyte        = c("Ly6c1", "Ccr2", "Ms4a6c", "F13a1", "Plbd1"),
  Macrophage      = c("Adgre1", "C1qa", "C1qb", "Cd68", "Siglec1", "Spic", "Timd4", "Mrc1"),
  DC              = c("Lamp3", "Clec4a3"),
  `B cell`        = c("Cxcr5"),
  Plasma          = c("Sdc1", "Slamf7"),
  `T cell`        = c("Cd3e", "Cd3d", "Cd4", "Cd8a", "Satb1", "Ctla4", "Foxp3", "Il2ra", "Pdcd1"),
  `NK/ILC`        = c("Ncr1", "Klrb1c", "Klra8"),
  MSC             = c("Col1a1", "Col1a2", "Pdgfra", "Pdgfrb", "Vcam1"),
  Endothelial     = c("Kdr", "Plvap", "Emb", "Tie1", "Clec14a", "Lyve1"),
  `Chondro/Osteo` = c("Col2a1", "Acan", "Comp", "Ucma", "Runx2")
)

# keep genes on the panel, each gene once (first group wins)
present = lapply(markers, intersect, rownames(mouse))
seen = character(0)
for (g in names(present)) { present[[g]] = setdiff(present[[g]], seen); seen = c(seen, present[[g]]) }
present = present[lengths(present) > 0]
setdiff(unlist(markers), rownames(mouse))   # markers not on your panel

ann_cols = c("celltype_seurat", "celltype_singleR",
  'celltype_seurat_cross', 'celltype_singleR_cross')
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

for (a in ann_cols) {
  p = DotPlot(mouse[,!is.na(mouse[[a]][,1])], features = present, group.by = a) + RotatedAxis()
  dir.create('../figures/bone_marrow/annotation', showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0('../figures/bone_marrow/annotation/mouse_dotplot_', sub('celltype_', '', a), '.png'), p,
         width = 24, height = 6, bg = 'white')
}

trusted = list(
celltype_seurat = c("MkP", "Granulocyte"),
celltype_singleR  = c("Monocyte", "MkP", "HSPC", "Granulocyte", 'Erythroid'),
celltype_seurat_cross = c('Plasma', 'MkP', 'Erythroid'),
celltype_singleR_cross = c('Plasma', 'MSC', 'MkP', 'Erythroid')
)
frac = 5   # chosen by marker AUROC (cell split): data/02_annotated_data/mouse_ebcc_frac_auroc.csv
n_k = list(
  celltype_seurat = table(mouse$celltype_seurat),
  celltype_singleR = table(mouse$celltype_singleR),
  celltype_seurat_cross = table(mouse$celltype_seurat_cross),
  celltype_singleR_cross = table(mouse$celltype_singleR_cross)
)
classes = sort(unique(na.omit(unlist(lapply(ann_cols, function(a) as.character(mouse@meta.data[[a]]))))))
make_beta = function(M) {
  beta = make_beta_prior(classes, ann_cols)
  # skip trusted classes the annotator never called (e.g. Plasma: not in the mouse reference)
  for (a in names(trusted)) for (k in intersect(trusted[[a]], names(n_k[[a]])[n_k[[a]] > 0]))
    beta[[a]][k, k] = 4 + frac * n_k[[a]][[k]] / M
  beta
}

# fit with M subtypes, then score how well the hidden labels are predicted
heldout_ll = function(M, R = 3) {
  f = run_ebcc(meta_ho, ann_cols, M = M, R = R, keep_rho = TRUE, beta = make_beta(M))
  lp = unlist(lapply(ann_cols, function(a) {
    s = ho[ho$ann == a, ]
    l = match(s$truth, f$classes); ok = !is.na(l)
    # P(hidden label) = sum over class-subtypes of P(cell in it) x P(annotator gives that label)
    log(rowSums(f$rho[s$cell[ok], , drop = FALSE] *
                t(f$par$v[[a]][, l[ok], drop = FALSE])))
  }))
  lp
}

Ms = c(1, 2, 3, 5, 10)
ll = sapply(Ms, heldout_ll)
mean(ll[,1])
mean(ll[,2])
mean(ll[,3])
mean(ll[,4])
mean(ll[,5])
mean(ll[,4]) - sd(ll[,4])/sqrt(dim(ll)[1])

# ---- fit ----
fit = run_ebcc(mouse@meta.data, ann_cols, M = 3, R = 10, beta = make_beta(3))

# ---- diagnostics ----
# 1. restarts: distance of each restart's ELBO from the best (large spread = many local optima)
round(fit$elbo_restarts - max(fit$elbo_restarts), 1)

# ---- save ----
mouse$ebcc_label   = fit$label
mouse$ebcc_maxpost = rowMaxs(fit$posterior)
mouse$ebcc_entropy = -rowSums(fit$posterior * log(fit$posterior + 1e-12))
# mouse$celltype   = mouse$ebcc_label   # switch over only once the diagnostics look sound

VariableFeatures(mouse) = rownames(mouse)
mouse = mouse |> 
NormalizeData()  |> 
ScaleData()
mouse$celltype = mouse$ebcc_label
saveRDS(fit$par,       'data/02_annotated_data/mouse_ebcc_params.rds')
saveRDS(fit$posterior, 'data/02_annotated_data/mouse_ebcc_posterior.rds')
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
as.anndata(mouse, file_path = 'data/02_annotated_data/',
           file_name = '02_mouse_annotated.h5ad', assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))

# marker dot plot for the final EBCC labels (mouse panel markers, `present` above)
dir.create('../figures/bone_marrow/annotation', showWarnings = FALSE, recursive = TRUE)
p = DotPlot(mouse, features = present, group.by = 'ebcc_label') + RotatedAxis() +
  labs(title = 'Mouse - EBCC labels', x = NULL, y = NULL)
ggsave('../figures/bone_marrow/annotation/mouse_dotplot_ebcc.png', p,
       width = 24, height = 6, bg = 'white')
