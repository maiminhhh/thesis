setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
human = readRDS('data/01_processed_data/01_human_processed.rds')
# references are keyed by Ensembl ID
rownames(human) = human@assays$Xenium@meta.data$ensembl_id

# background QC before annotation
human_cells = fread("data/00_raw_data/Xenium_V1_hColon_Non_diseased_Add_on_FFPE_outs/cells.csv.gz",
                    select = c("cell_id", "cell_area"))
human$area = human_cells$cell_area[match(colnames(human), human_cells$cell_id)]
keep_cells = bg_filter(nCount_probe_controls = human$nCount_ControlProbe,
                       nCount_transcripts = human$nCount_Xenium,
                       cell_areas = human$area,
                       n_genes = dim(human)[1], n_total_probes = 20)
human = human[,keep_cells]

human = human |>
  NormalizeData() |>
  ScaleData() |>
  RunPCA(features = rownames(human)) |>
  RunUMAP(dims = 1:20, return.model = TRUE) |>
  FindNeighbors(reduction = 'pca', dims = 1:20)
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done QC!')

# within-species reference (built in 02_human_annotate.qmd)
human_ref = readRDS('data/02_annotation_refs/human_processed_human_use.rds')

human_seurat_label = LabelTransfer(data = human, ref = human_ref)
human$celltype_seurat = human_seurat_label$cell_type_transfer
human$labelscore_seurat = human_seurat_label$score_transfer
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done Label Transfer!')

human_singleR_label = LabelSingleR(data = human, ref = human_ref)
human$celltype_singleR = human_singleR_label$cell_type_singleR
human$labeldelta_singleR = human_singleR_label$delta_singleR
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done SingleR!')

# cross-species reference
# ref source: https://www.parsebiosciences.com/datasets/5-million-mouse-single-cell-atlas-from-7-tissues/
mouse_ref_adata = read_h5ad('data/02_annotation_refs/mouse.h5ad')
mouse_count_matrix = t(mouse_ref_adata$X)
colnames(mouse_count_matrix) = mouse_ref_adata$obs$cell_ids
rownames(mouse_count_matrix) = mouse_ref_adata$var |> rownames()
mouse_metadata = as.data.frame(mouse_ref_adata$obs)
# the obs index (0, 1, ..., "0-1") is not cell_ids, and CreateSeuratObject matches
# meta.data rows to cells by name: name them like the count matrix columns
rownames(mouse_metadata) = mouse_ref_adata$obs$cell_ids
stopifnot(identical(rownames(mouse_metadata), colnames(mouse_count_matrix)))
mouse_metadata$barcode = mouse_metadata$barcode |>
  str_remove('_5d84e6de$')
mouse_ref = CreateSeuratObject(counts = mouse_count_matrix,
                               meta.data = mouse_metadata)
rm(mouse_ref_adata, mouse_count_matrix)
mouse_ref = mouse_ref |>
  NormalizeData()

homology_tbl = read.csv('../homology_tbl.csv')
# concatenate() already renames the mouse genes to their human orthologs;
# renaming again from the Xenium feature metadata scrambles the gene names
mouse_ref_human = concatenate(human_data = human,
                              mouse_data = mouse_ref, homology_tbl = homology_tbl,
                              tech_mouse = 'RNA')
DefaultAssay(mouse_ref_human) = 'RNA'
DefaultLayer(mouse_ref_human[['RNA']]) = 'data'
mouse_ref_human = mouse_ref_human |>
  subset(orig.ident == 'mmusculus')

mouse_ref_human$nCount_RNA = colSums(x = mouse_ref_human[['RNA']]$counts)
mouse_ref_human = mouse_ref_human |>
  subset(nCount_RNA > 0)

mouse_ref_human@meta.data = mouse_ref_human@meta.data |>
  mutate(across(c(B.Cell, Colonocyte.DCC, Enteric.Glial, Enteric.Neuron, Enterocyte,
                  Fibroblast, Goblet, Immune, PCC, Mesothelial, `Plasma.Cell`,
                  `Smooth.Muscle`, `Stem.TA`, Tuft, `Vascular.Endothelial`, Unknown),
                ~ na_if(.x, 'None'))) |>
  mutate(cell_type_detailed = coalesce(B.Cell, Colonocyte.DCC, Enteric.Glial,
                                       Enterocyte, Enteric.Neuron, Fibroblast, Goblet,
                                       Immune, PCC, Mesothelial, Plasma.Cell,
                                       Smooth.Muscle, Stem.TA, Tuft,
                                       Vascular.Endothelial, Unknown))
mouse_ref_human = mouse_ref_human |>
  subset(!is.na(cell_type_detailed) & !cell_type_detailed == 'Unknown')

mouse_cells = mouse_ref_human$cell_type_detailed
mouse_ref_human$cell_type_regroup = case_when(
  mouse_cells %in% c('Colonocyte/DCC', 'PCC', 'Enterocyte') ~ 'Epithelial: Absorptive',

  mouse_cells %in% c('Goblet', 'Tuft') ~ 'Epithelial: Secretory',

  mouse_cells %in% c('Stem/TA') ~ 'Epithelial: Stem/Progenitor',

  mouse_cells %in% c('Vascular Endothelial') ~ 'Endothelial',

  mouse_cells %in% c('Fibroblast', 'Mesothelial') ~ 'Stromal/Mesenchymal',

  mouse_cells %in% c('Smooth Muscle') ~ 'Smooth muscle',

  mouse_cells %in% c('Enteric Glial', 'Enteric Neuron') ~ 'Neural/Glial',

  mouse_cells %in% c('B Cell', 'Plasma Cell') ~ 'B/Plasma',

  mouse_cells %in% c('Immune') ~ 'Immune'
)

SaveSeuratRds(mouse_ref_human, 'data/02_annotation_refs/mouse_processed_human_use.rds')
print('Done pre-processing cross-species reference!')

human_ref = mouse_ref_human
human_seurat_label = LabelTransfer(data = human, ref = human_ref)
human$celltype_seurat_cross = human_seurat_label$cell_type_transfer
human$labelscore_seurat_cross = human_seurat_label$score_transfer
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done Label Transfer - Cross-species Reference!')

human_singleR_label = LabelSingleR(data = human, ref = human_ref)
human$celltype_singleR_cross = human_singleR_label$cell_type_singleR
human$labeldelta_singleR_cross = human_singleR_label$delta_singleR
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
print('Done SingleR - Cross-species Reference!')

# EBCC consensus
print('EBCC')
human = readRDS('data/02_annotated_data/02_human_annotated.rds')
rownames(human) = human[['Xenium']]@meta.data$gene_name

markers = list(
  `Epithelial: Absorptive` = c(
    "CA1", "CA2", "CA4", "SLC26A3", "GUCA2A", "GUCA2B",
    "MS4A12", "CEACAM7", "AQP8", "FABP1", "KRT20", "SELENBP1"
  ),
  `Epithelial: Secretory` = c(
    # Goblet
    "MUC2", "TFF3", "FCGBP", "SPINK4", "ZG16", "CLCA1", "AGR2", "REG4",
    # Enteroendocrine
    "CHGA", "CHGB", "NEUROD1", "PYY",
    # Tuft
    "POU2F3", "TRPM5", "DCLK1", "AVIL"
  ),
  `Epithelial: Stem/Progenitor` = c(
    "LGR5", "OLFM4", "ASCL2", "SMOC2", "RGMB", "SOX9",
    "MKI67", "TOP2A", "PCNA"
  ),
  Endothelial = c(
    "PECAM1", "CDH5", "VWF", "CLDN5", "PLVAP", "KDR",
    "FLT1", "EGFL7", "CLEC14A", "RAMP2"
  ),
  `Stromal/Mesenchymal` = c(
    "COL1A1", "COL1A2", "COL3A1", "COL6A2", "PDGFRA",
    "DCN", "LUM", "DPT", "FBLN1", "THY1"
  ),
  `Smooth muscle` = c(
    "ACTA2", "MYH11", "TAGLN", "DES", "MYLK",
    "CNN1", "ACTG2", "PLN", "DSTN"
  ),
  `Neural/Glial` = c(
    "S100B", "PLP1", "SOX10", "GFAP", "ELAVL4",
    "STMN2", "NRXN1", "KCNMB4"
  ),
  `B/Plasma` = c(
    # B cells
    "MS4A1", "CD19", "CD79A", "CD79B", "BANK1",
    # Plasma cells
    "MZB1", "JCHAIN", "DERL3", "XBP1", "PRDM1", "SDC1", "IGHA1", "IGHG1"
  ),
  Immune = c(
    "PTPRC",
    # T / NK
    "CD3D", "CD3E", "CD2", "CD8A", "NKG7",
    # Myeloid
    "LYZ", "CD68", "CD14", "C1QA", "ITGAM",
    # Mast
    "TPSAB1", "CPA3", "MS4A2"
  ),
  Erythroid = c(
    # membrane / cytoskeleton
    "SLC4A1", "SPTA1", "SPTB", "ANK1", "EPB41", "EPB42",
    "RHAG", "RHCE", "GYPB", "GYPC", "DMTN", "ADD2", "TMOD1",
    # heme / iron biosynthesis
    "ALAS2", "SLC25A37", "FECH", "UROD", "BLVRB", "AHSP",
    # transcription factors (less specific — don't use alone)
    "KLF1", "GATA1", "TAL1", "HEMGN"
  )
)

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

# chosen from the dot plots of the earlier run, when the cross-species reference was scrambled:
# re-check them against the new dot plots
trusted = list(
celltype_seurat = c('Epithelial: Absorptive', 'Endothelial', 'Stromal/Mesenchymal', 'Smooth muscle', 'B/Plasma'),
celltype_singleR  = c('Epithelial: Absorptive', 'Endothelial', 'Stromal/Mesenchymal', 'Smooth muscle', 'Neural/Glial', 'B/Plasma', 'Immune'),
celltype_seurat_cross = c('Epithelial: Absorptive', 'Endothelial', 'Stromal/Mesenchymal', 'Smooth muscle', 'B/Plasma'),
celltype_singleR_cross = c('Epithelial: Absorptive', 'Endothelial', 'Stromal/Mesenchymal', 'Smooth muscle', 'B/Plasma', 'Immune')
)
frac = 0.25   # chosen by marker AUROC with M fixed, 50/50 split (see human_ebcc_frac_auroc_chosen.csv)
n_k = lapply(setNames(ann_cols, ann_cols), function(a) table(human@meta.data[[a]]))
# every label any annotator gives, so the prior covers all classes run_ebcc sees
all_classes = human@meta.data[ann_cols] |> lapply(as.character) |> unlist() |> na.omit() |> unique()

# fit with M subtypes, then score how well the hidden labels are predicted
heldout_ll = function(M, R = 3) {
  beta = make_beta_prior(all_classes, ann_cols)
  for (a in names(trusted)) {for (k in trusted[[a]]){
  beta[[a]][k, k] = 4 + frac * n_k[[a]][[k]] / M }}
  f = run_ebcc(meta_ho, ann_cols, M = M, R = R, keep_rho = TRUE, beta = beta)
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
print(data.frame(M = Ms, heldout_ll = colMeans(ll),
                 se = apply(ll, 2, sd) / sqrt(nrow(ll))))

# ---- fit ----
M = 5   # chosen by the held-out (M, frac) grid, 1-SE rule (see human_ebcc_heldout_grid.csv)
beta = make_beta_prior(all_classes, ann_cols)
for (a in names(trusted)) for (k in trusted[[a]])
beta[[a]][k, k] = 4 + frac * n_k[[a]][[k]] / M
fit = run_ebcc(human@meta.data, ann_cols, M = M, R = 10, beta = beta)

# ---- diagnostics ----
# 1. restarts: distance of each restart's ELBO from the best (large spread = many local optima)
print(round(fit$elbo_restarts - max(fit$elbo_restarts), 1))

# ---- save ----
human$ebcc_label   = fit$label
human$ebcc_maxpost = rowMaxs(fit$posterior)
human$ebcc_entropy = -rowSums(fit$posterior * log(fit$posterior + 1e-12))
human$celltype = human$ebcc_label
print(table(human$celltype))

save_marker_dotplots(human, markers, c(ann_cols, 'ebcc_label'), prefix = 'human',
                     out_dir = '../figures/colon/annotation',
                     file_names = c('seurat', 'singleR', 'seurat_cross', 'singleR_cross', 'ebcc'))

VariableFeatures(human) = rownames(human)
human = human |>
NormalizeData()  |>
ScaleData()
saveRDS(fit$par,       'data/02_annotated_data/human_ebcc_params.rds')
saveRDS(fit$posterior, 'data/02_annotated_data/human_ebcc_posterior.rds')
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
as.anndata(human, file_path = 'data/02_annotated_data/',
           file_name = '02_human_annotated.h5ad', assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))
