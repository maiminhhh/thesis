setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'fetal_liver'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 16)   # SingleR
annotation_dir = '../figures/fetal_liver/annotation'   # dotplots and the EBCC grid heatmap
human = readRDS('data/01_processed_data/01_human_processed.rds')
VariableFeatures(human) = rownames(human)
human = human |>
    NormalizeData() |>
    ScaleData()
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done setting up!')

# within species reference
# ref source: https://developmental.cellatlas.io/fetal-liver
human_ref = zellkonverter::readH5AD('data/02_annotation_refs/human_ref.h5ad', reader = 'R')
human_ref = human_ref |>
    as.Seurat(counts = 'X', data = NULL)

human_ref[["RNA"]] = as(human_ref[["originalexp"]], "Assay5")
DefaultAssay(human_ref) = 'RNA'
human_ref = DietSeurat(human_ref, assays = 'RNA')

human_ref = human_ref |>
    NormalizeData()

human_cells = human_ref$cell.labels |> as.character()
human_ref$celltype = case_when(
    human_cells %in% c('HSC_MPP', 'MEMP') ~ 'HSPC',
    human_cells == 'Early Erythroid' ~ 'Erythroid progenitor',
    human_cells == 'Mid Erythroid' ~ 'Erythroblast',
    human_cells == 'Late Erythroid' ~ 'Erythrocyte',
    human_cells == 'Megakaryocyte' ~ 'Megakaryocyte',
    human_cells %in% c('Monocyte', 'Kupffer Cell', 'Mono-Mac', 'Monocyte precursor', 'VCAM1+ EI macrophage') ~ 'Mono/Macrophage',
    human_cells == 'Neutrophil-myeloid progenitor' ~ 'Granulocyte',
    human_cells == 'Hepatocyte' ~ 'Hepatocyte',
    human_cells == 'Endothelial cell' ~ 'Endothelial',
    human_cells == 'Mast cell' ~ 'Mast cell',
    human_cells %in% c('DC1', 'DC2', 'DC precursor', 'pDC precursor') ~ 'Dendritic cell',
    human_cells %in% c('Pre pro B cell', 'pro-B cell', 'pre-B cell', 'B cell') ~ 'B cell',
    human_cells %in% c('Early lymphoid_T lymphocyte', 'ILC precursor', 'NK') ~ 'T/NK/ILC',
    human_cells == 'Fibroblast' ~ 'Fibroblast'
) |>  as.factor()
stopifnot(!anyNA(human_ref$celltype))
SaveSeuratRds(human_ref, 'data/02_annotation_refs/02_human_ref.rds')
print('Done pre-processing reference!')

human_intersect = intersect(rownames(human), rownames(human_ref))
human_annot = human[human_intersect,]
human_ref = human_ref[human_intersect,]

# label transfer
human_seurat_label = LabelTransfer(data = human_annot, ref = human_ref,
                                   celltype_column = 'celltype')
human$celltype_seurat = human_seurat_label$cell_type_transfer |> as.factor()
human$labelscore_seurat = human_seurat_label$score_transfer
print('Done Label Transfer!')

# singleR
human_singleR_label = LabelSingleR(data = human_annot, ref = human_ref,
                                   celltype_column = 'celltype',
                                   tech = 'Xenium', BPPARAM = param)
human$celltype_singleR = human_singleR_label$cell_type_singleR |>
    as.factor()
human$labeldelta_singleR = human_singleR_label$delta_singleR
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done SingleR!')

# cross-species reference
mouse_ref = readRDS('data/02_annotation_refs/02_mouse_ref.rds')
homology_tbl = read.csv('../homology_tbl.csv')
mouse_ref_human = cross_species_ref(mouse_ref, human, homology_tbl,
                                    ref_species = 'mmusculus')
mouse_ref_human = mouse_ref_human |>
    subset(nCount_RNA > 0)
SaveSeuratRds(mouse_ref_human, 'data/02_annotation_refs/mouse_processed_human_use.rds')
print('Done pre-processing cross-species reference!')

human_intersect = intersect(rownames(human), rownames(mouse_ref_human))
human_annot = human[human_intersect,]
human_ref = mouse_ref_human[human_intersect,]

# label transfer
human_seurat_label = LabelTransfer(data = human_annot, ref = human_ref,
                                   celltype_column = 'celltype')
human$celltype_seurat_cross = human_seurat_label$cell_type_transfer |> as.factor()
human$labelscore_seurat_cross = human_seurat_label$score_transfer
print('Done Label Transfer - Cross-species Reference!')

# singleR
human_singleR_label = LabelSingleR(data = human_annot, ref = human_ref,
                                   celltype_column = 'celltype',
                                   tech = 'Xenium', BPPARAM = param)
human$celltype_singleR_cross = human_singleR_label$cell_type_singleR |>
    as.factor()
human$labeldelta_singleR_cross = human_singleR_label$delta_singleR
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done SingleR - Cross-species Reference!')

# marker dot plots, one per annotator
# every gene is on the human Xenium panel (the canonical ALB/HB*/PF4-type markers are not)
markers = list(
  "Hepatocyte"           = c("GPC3", "APOE", "HAMP", "TAT", "HPX", "CRP", "PLG", "KNG1", "SELENOP",
                             "HMGCS2", "ADH4", "ADH1C", "APOA5", "CFHR1", "CYP3A4"),
  "HSPC"                 = c("CD34", "PROM1", "CRHBP", "HOXA9"),
  "Erythroid progenitor" = c("GATA1", "KLF1", "TFRC", "CA1", "EPOR"),
  "Erythroblast"         = c("GYPA", "GYPB", "ALAS2", "SLC4A1", "AHSP", "HEMGN", "RHAG",
                             "MKI67", "TOP2A"),
  "Erythrocyte"          = c("SNCA", "CA2", "ICAM4", "BNIP3L"),
  "Megakaryocyte"        = c("GNG11", "CAVIN2", "VWF", "MMRN1"),
  "Mast cell"            = c("CPA3", "MS4A2", "KIT", "GATA2", "HPGDS", "IL1RL1", "SLC18A2", "VWA5A"),
  "Granulocyte"          = c("MPO", "ELANE", "CTSG", "CAMP", "S100A12", "MNDA"),
  "Mono/Macrophage"      = c("CD14", "LYZ", "FCN1", "VCAN", "FCGR3A", "CD68", "AIF1", "CD163",
                             "C1QA", "C1QB", "C1QC", "MARCO", "VCAM1", "CD5L", "FOLR2", "LYVE1",
                             "MRC1", "SIGLEC1"),
  "Dendritic cell"       = c("IRF8", "CD1C", "CD1E", "CD1A", "FCER1A", "CLEC10A", "LILRA4", "IL3RA",
                             "PLD4", "LAMP3", "CCR7", "FLT3"),
  "B cell"               = c("CD79A", "CD19", "PAX5", "MS4A1", "BANK1"),
  "T/NK/ILC"             = c("CD3E", "CD3D", "CD2", "CD247", "CD8A", "IL7R", "KLRB1", "NKG7",
                             "GNLY", "KLRD1", "GZMA", "GZMB", "PRF1"),
  "Endothelial"          = c("PECAM1", "KDR", "EGFL7", "CLEC14A", "RAMP2", "ERG", "FCN2",
                             "DNASE1L3", "ACKR1"),
  "Fibroblast"           = c("COL1A1", "COL1A2", "BGN", "PDGFRA", "PDGFRB", "PCOLCE", "OGN", "DPT",
                             "FBLN1", "SFRP2", "RBP5", "DES", "ACTA2", "MYH11")
)

ann_cols = c("celltype_seurat", "celltype_singleR",
  'celltype_seurat_cross', 'celltype_singleR_cross')
save_marker_dotplots(human, markers, ann_cols, prefix = '02_human', out_dir = annotation_dir)
print('Done dot plots!')

print('EBCC')
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
celltype_seurat = c("Hepatocyte", 'Megakaryocyte', "Mast cell", "Granulocyte", 'Mono/Macrophage', 'Dendritic cell', 'B cell', 'T/NK/ILC', 'Endothelial', 'Fibroblast'),
celltype_singleR  = c("Hepatocyte", 'Megakaryocyte', "Mast cell", "Granulocyte", 'Mono/Macrophage', 'B cell', 'T/NK/ILC', 'Endothelial', 'Fibroblast'),
celltype_seurat_cross = c('Mono/Macrophage'),
celltype_singleR_cross = c('Endothelial', 'Hepatocyte')
)

n_k = list(
  celltype_seurat = table(human$celltype_seurat),
  celltype_singleR = table(human$celltype_singleR),
  celltype_seurat_cross = table(human$celltype_seurat_cross),
  celltype_singleR_cross = table(human$celltype_singleR_cross)
)
# every label any annotator gives, so no annotator's label falls outside the class set
classes = sort(unique(na.omit(unlist(lapply(ann_cols, function(a) as.character(human@meta.data[[a]]))))))
make_beta = function(M, frac) {
  beta = make_beta_prior(classes, ann_cols)
  # skip trusted classes the annotator never called
  for (a in names(trusted)) for (k in intersect(trusted[[a]], names(n_k[[a]])[n_k[[a]] > 0]))
    beta[[a]][k, k] = 4 + frac * n_k[[a]][[k]] / M
  beta
}

# choose M and frac (ebcc_grid.R, ../EBCC_PIPELINE.md sections 3 and 3b)
source('ebcc_grid.R')
Ms = c(1, 2, 3, 5, 10)
fracs = c(0, 0.25, 0.5, 1, 2)
# last chosen values, to report a change: data/02_annotated_data/human_ebcc_frac_auroc_chosen.csv
ebcc_previous = list(M = 5, frac = 0)
# M: held-out log-likelihood over the (M, frac) grid; among (M, frac) within 1 SE of the
# best, the smallest M. Only M is taken from this grid
run_ebcc_grid('human', Ms, fracs, meta_ho, ho, ann_cols, make_beta, trusted)
M = select_ebcc_grid('human', Ms, fracs, ebcc_previous, fig_dir = annotation_dir)$M
# frac: with M fixed, a full fit (R = 10) for each frac; the frac with the highest mean
# marker AUROC on a random half of the cells (reported on the other half). frac is tuned
# for marker agreement, so the marker AUROC is not an independent check of EBCC
fracs_auroc = c(0, 0.25, 0.5, 1, 2, 5)
run_ebcc_frac_fits('human', human@meta.data, ann_cols, M, fracs_auroc, make_beta, trusted)
choice = select_frac_by_auroc('human', human, M, fracs_auroc, markers, ebcc_previous,
                              fig_dir = annotation_dir)
frac = choice$frac
# EBCC_FIT_IF_CHANGED=1: keep the existing labels when neither M nor frac changed
if (Sys.getenv('EBCC_FIT_IF_CHANGED') == '1' && !choice$changed) quit(save = 'no')

# ---- fit: the saved full fit for the chosen frac ----
fit = choice$fit

# ---- diagnostics ----
# 1. restarts: distance of each restart's ELBO from the best (large spread = many local optima)
print(round(fit$elbo_restarts - max(fit$elbo_restarts), 1))

# ---- save ----
human$ebcc_label   = fit$label
human$ebcc_maxpost = rowMaxs(fit$posterior)
human$ebcc_entropy = -rowSums(fit$posterior * log(fit$posterior + 1e-12))
human$celltype = human$ebcc_label
print(table(human$celltype))

save_marker_dotplots(human, markers, 'ebcc_label', prefix = 'human', out_dir = annotation_dir)
file.rename(file.path(annotation_dir, 'human_dotplot_ebcc_label.png'), file.path(annotation_dir, 'human_dotplot_ebcc.png'))

saveRDS(fit$par,       'data/02_annotated_data/human_ebcc_params.rds')
saveRDS(fit$posterior, 'data/02_annotated_data/human_ebcc_posterior.rds')
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
as.anndata(human, file_path = 'data/02_annotated_data/',
           file_name = '02_human_annotated.h5ad', assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))
