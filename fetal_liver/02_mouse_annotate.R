setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'fetal_liver'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
param = MulticoreParam(workers = 16)   # SingleR
annotation_dir = '../figures/fetal_liver/annotation'   # dotplots and the EBCC grid heatmap
mouse = readRDS('data/01_processed_data/01_mouse_processed.rds')
VariableFeatures(mouse) = rownames(mouse)
mouse = mouse |>
    NormalizeData() |>
    ScaleData()
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
print('Done setting up!')

# within species reference
# ref source: https://singlocell.openscience.mcgill.ca/display?dataset=RNA_Ms_Liver_2021
mouse_ref = readRDS('data/02_annotation_refs/GSE180050_combined.rds')
# default assay is 'integrated' (2000 HVGs, no counts): use the full RNA assay
DefaultAssay(mouse_ref) = 'RNA'
mouse_ref = DietSeurat(mouse_ref, assays = 'RNA')
mouse_ref[['RNA']] = as(mouse_ref[['RNA']], 'Assay5')
mouse_ref = mouse_ref |>
    NormalizeData()
mouse_cells = mouse_ref$CellType |>
    as.character()
mouse_ref$celltype = case_when(
    mouse_cells == 'HPSC' ~ 'HSPC',
    mouse_cells == 'EryP' ~ 'Erythroid progenitor',
    mouse_cells == 'Eryb' ~ 'Erythroblast',
    mouse_cells == 'Ery' ~ 'Erythrocyte',
    mouse_cells == 'Mk' ~ 'Megakaryocyte',
    mouse_cells == 'Mo/KC' ~ 'Mono/Macrophage',
    mouse_cells == 'G' ~ 'Granulocyte',
    mouse_cells == 'Hepa' ~ 'Hepatocyte',
    mouse_cells == 'EC' ~ 'Endothelial'
)  |> factor()
stopifnot(!anyNA(mouse_ref$celltype))
SaveSeuratRds(mouse_ref, 'data/02_annotation_refs/02_mouse_ref.rds')
print('Done pre-processing reference!')

mouse_intersect = intersect(rownames(mouse), rownames(mouse_ref))
mouse_ref = mouse_ref[mouse_intersect,]
mouse_annot = mouse[mouse_intersect,]

# label transfer
mouse_seurat_label = LabelTransfer(data = mouse_annot, ref = mouse_ref,
                                   celltype_column = 'celltype')
mouse$celltype_seurat = mouse_seurat_label$cell_type_transfer |> as.factor()
mouse$labelscore_seurat = mouse_seurat_label$score_transfer
print('Done Label Transfer!')

# singleR
mouse_singleR_label = LabelSingleR(data = mouse_annot, ref = mouse_ref,
                                   celltype_column = 'celltype',
                                   tech = 'Xenium', BPPARAM = param)
mouse$celltype_singleR = mouse_singleR_label$cell_type_singleR |>
    as.factor()
mouse$labeldelta_singleR = mouse_singleR_label$delta_singleR
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
print('Done SingleR!')

# cross-species reference
human_ref = readRDS('data/02_annotation_refs/02_human_ref.rds')
homology_tbl = read.csv('../homology_tbl.csv')
human_ref_mouse = cross_species_ref(human_ref, mouse, homology_tbl,
                                    ref_species = 'hsapiens')
human_ref_mouse = human_ref_mouse |>
    subset(nCount_RNA > 0)
SaveSeuratRds(human_ref_mouse, 'data/02_annotation_refs/human_processed_mouse_use.rds')
print('Done pre-processing cross-species reference!')

mouse_intersect = intersect(rownames(mouse), rownames(human_ref_mouse))
mouse_annot = mouse[mouse_intersect,]
mouse_ref = human_ref_mouse[mouse_intersect,]

# label transfer
mouse_seurat_label = LabelTransfer(data = mouse_annot, ref = mouse_ref,
                                   celltype_column = 'celltype')
mouse$celltype_seurat_cross = mouse_seurat_label$cell_type_transfer |> as.factor()
mouse$labelscore_seurat_cross = mouse_seurat_label$score_transfer
print('Done Label Transfer - Cross-species Reference!')

# singleR
mouse_singleR_label = LabelSingleR(data = mouse_annot, ref = mouse_ref,
                                   celltype_column = 'celltype',
                                   tech = 'Xenium', BPPARAM = param)
mouse$celltype_singleR_cross = mouse_singleR_label$cell_type_singleR |>
    as.factor()
mouse$labeldelta_singleR_cross = mouse_singleR_label$delta_singleR
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
print('Done SingleR - Cross-species Reference!')


# marker dot plots, one per annotator
markers = list(
  "Hepatocyte"           = c("Ambp","Hpx","Trf","Gc","Igfbp1"),
  "HSPC"                 = c("Kit","Ly6a","Cd34","Slamf1","Hemgn"),
  "Erythroid progenitor" = c("Gata1","Klf1","Epor","Car2"),
  "Erythroblast"         = c("Gypa","Tfrc","Slc4a1","Rhd","Fech"),
  "Erythrocyte"          = c("Hbb-bt"),
  "Megakaryocyte"        = c("Pf4","Ppbp","Gng11"),
  "Mast cell"            = c("Cpa3","Ltc4s"),
  "Granulocyte"          = c("Mpo","Elane","Ctsg","Ly6g","S100a9"),
  "Mono/Macrophage"      = c("Adgre1","C1qa","Cd68","Clec4f","Vcam1"),
  "Dendritic cell"       = c("Plbd1","Lamp3","Clec4a3"),
  "B cell"               = c("Cd93","Cxcr5"),
  "T/NK/ILC"             = c("Cd3e","Cd3d","Ncr1","Klrb1c"),
  "Endothelial"          = c("Kdr","Clec14a","Tie1","Dnase1l3","Oit3"),
  "Fibroblast"           = c("Col1a1","Col1a2","Pdgfrb","Des","Colec11")
)

ann_cols = c("celltype_seurat", "celltype_singleR",
  'celltype_seurat_cross', 'celltype_singleR_cross')
save_marker_dotplots(mouse, markers, ann_cols, prefix = '02_mouse', out_dir = annotation_dir)
print('Done dot plots!')

print('EBCC')
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

trusted = list(
celltype_seurat = c("Hepatocyte", "Granulocyte", 'Mono/Macrophage', 'Endothelial'),
celltype_singleR  = c("Hepatocyte", "HSPC", "Erythroblast", 'Megakaryocyte', 'Granulocyte', 'Mono/Macrophage', 'Endothelial'),
celltype_seurat_cross = c('Hepatocyte', 'Megakaryocyte', 'Mast cell', 'Granulocyte', 'Mono/Macrophage', 'Endothelial', 'Fibroblast'),
celltype_singleR_cross = c('Hepatocyte', 'HSPC', 'Mast cell', 'Granulocyte', 'Mono/Macrophage', 'Endothelial', 'Fibroblast')
)
n_k = list(
  celltype_seurat = table(mouse$celltype_seurat),
  celltype_singleR = table(mouse$celltype_singleR),
  celltype_seurat_cross = table(mouse$celltype_seurat_cross),
  celltype_singleR_cross = table(mouse$celltype_singleR_cross)
)
classes = sort(unique(na.omit(unlist(lapply(ann_cols, function(a) as.character(mouse@meta.data[[a]]))))))

make_beta = function(M, frac) {
  beta = make_beta_prior(classes, ann_cols)
  # skip trusted classes the annotator never called (e.g. Plasma: not in the mouse reference)
  for (a in names(trusted)) for (k in intersect(trusted[[a]], names(n_k[[a]])[n_k[[a]] > 0]))
    beta[[a]][k, k] = 4 + frac * n_k[[a]][[k]] / M
  beta
}

# choose M and frac (ebcc_grid.R, ../EBCC_PIPELINE.md sections 3 and 3b)
source('ebcc_grid.R')
Ms = c(1, 2, 3, 5, 10)
fracs = c(0, 0.25, 0.5, 1, 2)
# last chosen values, to report a change: data/02_annotated_data/mouse_ebcc_frac_auroc_chosen.csv
ebcc_previous = list(M = 3, frac = 0)
# M: held-out log-likelihood over the (M, frac) grid; among (M, frac) within 1 SE of the
# best, the smallest M. Only M is taken from this grid
run_ebcc_grid('mouse', Ms, fracs, meta_ho, ho, ann_cols, make_beta, trusted)
M = select_ebcc_grid('mouse', Ms, fracs, ebcc_previous, fig_dir = annotation_dir)$M
# frac: with M fixed, a full fit (R = 10) for each frac; the frac with the highest mean
# marker AUROC on a random half of the cells (reported on the other half). frac is tuned
# for marker agreement, so the marker AUROC is not an independent check of EBCC
fracs_auroc = c(0, 0.25, 0.5, 1, 2, 5)
run_ebcc_frac_fits('mouse', mouse@meta.data, ann_cols, M, fracs_auroc, make_beta, trusted)
choice = select_frac_by_auroc('mouse', mouse, M, fracs_auroc, markers, ebcc_previous,
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
mouse$ebcc_label   = fit$label
mouse$ebcc_maxpost = rowMaxs(fit$posterior)
mouse$ebcc_entropy = -rowSums(fit$posterior * log(fit$posterior + 1e-12))
mouse$celltype = mouse$ebcc_label
print(table(mouse$celltype))

save_marker_dotplots(mouse, markers, 'ebcc_label', prefix = 'mouse', out_dir = annotation_dir)
file.rename(file.path(annotation_dir, 'mouse_dotplot_ebcc_label.png'), file.path(annotation_dir, 'mouse_dotplot_ebcc.png'))

saveRDS(fit$par,       'data/02_annotated_data/mouse_ebcc_params.rds')
saveRDS(fit$posterior, 'data/02_annotated_data/mouse_ebcc_posterior.rds')
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
as.anndata(mouse, file_path = 'data/02_annotated_data/',
           file_name = '02_mouse_annotated.h5ad', assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))
