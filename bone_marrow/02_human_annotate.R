source('../library.R')
source('../source.R')
human = readRDS('data/01_processed_data/01_human_processed.rds')
human = human |> 
    NormalizeData() |> 
    ScaleData()
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done setting up!')

# ref source: https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79/project-matrices
human_ref = read_h5ad('data/02_annotation_refs/511478a9-1940-4aa3-ab63-e2791aa6e623/BM_pooling_and_control.h5ad')
human_ref = human_ref$as_Seurat()
# X is log1p(CP100k), not counts; LogNormalize is per-cell scale invariant,
# so expm1(X) normalizes identically to the raw counts
human_ref[['RNA']]$counts = expm1(human_ref[['RNA']]$X)
human_ref[['RNA']]$X = NULL

human_ref = human_ref |>
    NormalizeData()

human_celltype = human_ref$anno
human_ref$celltype = case_when(
    human_celltype == 'HSCs' ~ 'HSPC',
    human_celltype == 'Megakaryocyte progenitors' ~ 'MkP',
    human_celltype %in% c('Erythroid progenitors', 'Erythroid cells', 'ANK1-low erythroid cells') ~ 'Erythroid',
    human_celltype == 'Neutrophil progenitors' ~ 'Granulocyte',
    human_celltype %in% c('CD14+ monocytes', 'CD16+ monocytes') ~ 'Monocyte',
    human_celltype %in% c('cDCs', 'pDCs') ~ 'DC',
    human_celltype %in% c('Pro-B cells', 'Pre-B cells', 'B cells') ~ 'B cell',
    human_celltype %in% c('CD4+ naive T cells', 'CD8+ naive T cells', 'T helper cells', 'Cytotoxic T cells') ~ 'T cell',
    human_celltype == 'NK cells' ~ 'NK/ILC',
    human_celltype == 'Plasma cells' ~ 'Plasma',
    human_celltype == 'MSCs' ~ 'MSC'
)  |> factor()

SaveSeuratRds(human_ref, 'data/02_annotation_refs/02_human_ref.rds')
human_ref = readRDS('data/02_annotation_refs/02_human_ref.rds')
human = readRDS('data/02_annotated_data/02_human_annotated.rds')
rownames(human) = human[['Xenium']]@meta.data$gene_name
human_intersect = intersect(rownames(human), rownames(human_ref))
human_annot = human[human_intersect,]
human_ref = human_ref[human_intersect,]
print('Done pre-processing reference!')

# label transfer
human_seurat_label = LabelTransfer(data = human_annot, ref = human_ref,
                                   celltype_column = 'celltype')
human$celltype_seurat = human_seurat_label$cell_type_transfer |> as.factor()
human$labelscore_seurat = human_seurat_label$score_transfer
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done Label Transfer!')

# singleR
human_singleR_label = LabelSingleR(data = human_annot, ref = human_ref, 
                                   celltype_column = 'celltype',
                                   tech = 'Xenium') 
human$celltype_singleR = human_singleR_label$cell_type_singleR |>
    as.factor()
human$labeldelta_singleR = human_singleR_label$delta_singleR
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done SingleR!')

# cross-species reference
human = readRDS('data/02_annotated_data/02_human_annotated.rds')
mouse_ref = readRDS('data/02_annotation_refs/02_mouse_ref.rds')
rownames(human) = human[['Xenium']]@meta.data$gene_id
homology_tbl = read.csv('../homology_tbl.csv')
mouse_gene_dict = homology_tbl |> 
  filter(!is.na(mmusculus_homolog_orthology_type)) |> 
  dplyr::select(mmusculus_homolog_ensembl_gene,
                mmusculus_homolog_associated_gene_name) |> 
  filter(!is.na(mmusculus_homolog_associated_gene_name)) |> 
  filter(!mmusculus_homolog_associated_gene_name == '') |> 
  distinct()
mouse_ref[['RNA']]@meta.data$gene_name = rownames(mouse_ref)
mouse_ref[['RNA']]@meta.data$ens = mouse_gene_dict$mmusculus_homolog_ensembl_gene[
  match(mouse_ref[['RNA']]@meta.data$gene_name, 
  mouse_gene_dict$mmusculus_homolog_associated_gene_name)
]
mouse_ref[['RNA']]@meta.data$gene_id = coalesce(mouse_ref[['RNA']]@meta.data$ens,
                                                mouse_ref[['RNA']]@meta.data$gene_name)
rownames(mouse_ref) = mouse_ref[['RNA']]@meta.data$gene_id
mouse_ref_human = concatenate(human, mouse_ref, homology_tbl, tech_mouse = 'RNA')
DefaultAssay(mouse_ref_human) = 'RNA'
DefaultLayer(mouse_ref_human[['RNA']]) = 'data'
mouse_ref_human = mouse_ref_human |> 
    subset(orig.ident == 'mmusculus')

mouse_ref_human$nCount_RNA = colSums(x = mouse_ref_human[['RNA']]$counts)
mouse_ref_human = mouse_ref_human |>
    subset(nCount_RNA > 0)
SaveSeuratRds(mouse_ref_human, 'data/02_annotation_refs/mouse_processed_human_use.rds')
print('Done pre-processing cross-species reference!')

# label transfer
mouse_ref_human = readRDS('data/02_annotation_refs/mouse_processed_human_use.rds')
rownames(human) = human[['Xenium']]@meta.data$gene_id
human_intersect = intersect(rownames(human), rownames(mouse_ref_human))
human_annot = human[human_intersect,]
xref = mouse_ref_human[human_intersect,]

human_seurat_label = LabelTransfer(data = human_annot, ref = xref,
                                   celltype_column = 'celltype')
human$celltype_seurat_cross = human_seurat_label$cell_type_transfer |> as.factor()
human$labelscore_seurat_cross = human_seurat_label$score_transfer
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done Label Transfer - Cross-species Reference!')

# singleR
human_singleR_label = LabelSingleR(data = human_annot, ref = xref,
                                   celltype_column = 'celltype',
                                   tech = 'Xenium')
human$celltype_singleR_cross = human_singleR_label$cell_type_singleR |>
    as.factor()
human$labeldelta_singleR_cross = human_singleR_label$delta_singleR
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done SingleR - Cross-species Reference!')

print('EBCC')
human = readRDS('data/02_annotated_data/02_human_annotated.rds')
rownames(human) = human[['Xenium']]@meta.data$gene_name

markers = list(
  HSPC          = c("CD34", "SPINK2", "PROM1", "CRHBP"),
  MkP           = c("ITGA2B", "PF4", "GP9", "PPBP"),
  Erythroid     = c("HBB", "GYPA", "CA1", "AHSP", "KLF1"),
  Granulocyte   = c("MPO", "ELANE", "AZU1", "LTF", "CAMP", "CSF3R", "FCGR3B", "CXCR2"),
  Monocyte      = c("CD14", "FCN1", "VCAN", "CSF1R", "FCGR3A"),
  DC            = c("CLEC4C", "LILRA4", "IL3RA", "CD1C", "FCER1A", "CLEC9A"),
  `B cell`      = c("MS4A1", "CD79A", "CD19", "PAX5", "VPREB1"),
  Plasma        = c("MZB1", "JCHAIN", "SDC1", "TNFRSF17", "PRDM1"),
  `T cell`      = c("CD3E", "CD3D", "CD4", "CD8A", "IL7R"),
  `NK/ILC`      = c("NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1"),
  MSC           = c("CXCL12", "LEPR", "KITLG", "PDGFRA", "COL1A1")
)

# keep genes on the panel, each gene once (first group wins)
present = lapply(markers, intersect, rownames(human))
seen = character(0)
for (g in names(present)) { present[[g]] = setdiff(present[[g]], seen); seen = c(seen, present[[g]]) }
present = present[lengths(present) > 0]
setdiff(unlist(markers), rownames(human))   # markers not on your panel

ann_cols = c("celltype_seurat", "celltype_singleR", 
  'celltype_seurat_cross', 'celltype_singleR_cross')
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

DotPlot(human, features = markers, group.by = 'celltype_seurat') + RotatedAxis()
DotPlot(human[,!is.na(human$celltype_singleR)], features = markers, group.by = 'celltype_singleR') + RotatedAxis()
DotPlot(human, features = markers, group.by = 'celltype_seurat_cross') + RotatedAxis()
DotPlot(human[,!is.na(human$celltype_singleR_cross)], features = markers, group.by = 'celltype_singleR_cross') + RotatedAxis()
for (a in ann_cols) {
  p = DotPlot(human[,!is.na(human[[a]][,1])], features = markers, group.by = a) + RotatedAxis()
  dir.create('../figures/bone_marrow/annotation', showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0('../figures/bone_marrow/annotation/human_dotplot_', sub('celltype_', '', a), '.png'), p,
         width = 18, height = 6, bg = 'white')
}

trusted = list(
celltype_seurat = c("NK/ILC", 'Plasma', "MSC", "HSPC", 'Granulocyte', 'Erythroid', 'B cell', 'T cell'),
celltype_singleR  = c('HSPC', "T cell", "Plasma", "NK/ILC", "MSC", 'Monocyte', 'Granulocyte', 'Erythroid', 'B cell'),
celltype_seurat_cross = c('Erythroid', 'T cell'),
celltype_singleR_cross = c('Erythroid', 'T cell')
)
frac = 5   # chosen by marker AUROC (cell split): data/02_annotated_data/human_ebcc_frac_auroc.csv
n_k = list(
  celltype_seurat = table(human$celltype_seurat),
  celltype_singleR = table(human$celltype_singleR),
  celltype_seurat_cross = table(human$celltype_seurat_cross),
  celltype_singleR_cross = table(human$celltype_singleR_cross)
)
# every label any annotator gives, so no annotator's label falls outside the class set
classes = sort(unique(na.omit(unlist(lapply(ann_cols, function(a) as.character(human@meta.data[[a]]))))))
make_beta = function(M) {
  beta = make_beta_prior(classes, ann_cols)
  # skip trusted classes the annotator never called
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
mean(ll[,5]) - sd(ll[,5])/sqrt(dim(ll)[1])

# ---- fit ----
M = 5
fit = run_ebcc(human@meta.data, ann_cols, M = M, R = 10, beta = make_beta(M))

# ---- diagnostics ----
# 1. restarts: distance of each restart's ELBO from the best (large spread = many local optima)
round(fit$elbo_restarts - max(fit$elbo_restarts), 1)

# ---- save ----
human$ebcc_label  = fit$label
human$ebcc_maxpost = rowMaxs(fit$posterior)
human$ebcc_entropy = -rowSums(fit$posterior * log(fit$posterior + 1e-12))
# human$celltype   = human$ebcc_label   # switch over only once the diagnostics look sound

VariableFeatures(human) = rownames(human)
human = human |> 
NormalizeData()  |> 
ScaleData()
human$celltype = human$ebcc_label
saveRDS(fit$par,       'data/02_annotated_data/human_ebcc_params.rds')
saveRDS(fit$posterior, 'data/02_annotated_data/human_ebcc_posterior.rds')
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
as.anndata(human, file_path = 'data/02_annotated_data/',
           file_name = '02_human_annotated.h5ad', assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))

# marker dot plot for the final EBCC labels (panel markers, each gene once)
present = lapply(markers, intersect, rownames(human))
seen = character(0)
for (g in names(present)) { present[[g]] = setdiff(present[[g]], seen); seen = c(seen, present[[g]]) }
present = present[lengths(present) > 0]
dir.create('../figures/bone_marrow/annotation', showWarnings = FALSE, recursive = TRUE)
p = DotPlot(human, features = present, group.by = 'ebcc_label') + RotatedAxis() +
  labs(title = 'Human - EBCC labels', x = NULL, y = NULL)
ggsave('../figures/bone_marrow/annotation/human_dotplot_ebcc.png', p,
       width = 18, height = 6, bg = 'white')
