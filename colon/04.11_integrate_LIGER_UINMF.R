print('1. LIGER UINMF Set up')
print('--------------------------------------')
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
source('integration_datasets.R')
options(future.globals.maxSize = 20 * 1024^3)

for (d in names(datasets)) {
print(str_glue('===== {d} ====='))
data = readRDS(str_glue('{datasets[[d]]$concat}.rds'))
out_path = datasets[[d]]$out
# per-species inputs from 03: the cells kept in the concatenated data, all genes
human = readRDS(str_glue('{datasets[[d]]$inputs}02_human_annotated.rds'))
human$orig.ident = 'hsapiens'
mouse = readRDS(str_glue('{datasets[[d]]$inputs}02_mouse_annotated.rds'))
mouse$orig.ident = 'mmusculus'

print('2. LIGER UINMF Preprocessing')
print('--------------------------------------')
# shared genes take their human id; mouse-only genes keep their mouse symbol
orthologs_mapping = data[['Xenium']]@meta.data
mm_to_hs = setNames(orthologs_mapping$hsapiens_homolog_ensembl_gene,
                    orthologs_mapping$mmusculus_homolog_ensembl_gene)
mouse_genes = mouse[['Xenium']]@meta.data
rownames(mouse) = coalesce(unname(mm_to_hs[mouse_genes$ensembl_id]), mouse_genes$gene_name)
rownames(human) = human[['Xenium']]@meta.data$ensembl_id |>
  coalesce(human[['Xenium']]@meta.data$gene_name)
print(sprintf('shared genes: %d | human-only: %d | mouse-only: %d',
              length(intersect(rownames(human), rownames(mouse))),
              length(setdiff(rownames(human), rownames(mouse))),
              length(setdiff(rownames(mouse), rownames(human)))))

print('3. Run LIGER UINMF')
print('--------------------------------------')
# prefixed cell names ('hsapiens_<barcode>') so H rows match the concatenated data
data_liger_uinmf = integrate_liger_uinmf(human, mouse, add_prefix = 'auto')
stopifnot(all(grepl('^(hsapiens|mmusculus)_', colnames(data_liger_uinmf$seurat))))

print('4. LIGER UINMF Save data')
print('--------------------------------------')
saveRDS(data_liger_uinmf$liger, str_glue('{out_path}04.11_LIGER_UINMF.rds'))

# H and W first: 04.12-04.17 need them even if the h5ad export below fails
print('5. LIGER UINMF Extract H matrix')
print('--------------------------------------')
H = getMatrix(data_liger_uinmf$liger, 'H')
# cells x factors, each cell scaled to unit length. Cells loading on < 2 factors
# (very low counts) are dropped: RPCA UINMF returns NA for them, and every UINMF
# method and the evaluation (eval_cells.txt) must use the same cells
H = lapply(H, function(M) l2(t(M)))
H = lapply(H, function(h) h[rowSums(h != 0) >= 2, , drop = FALSE])
H_raw = rbind(H$mmusculus,
          H$hsapiens)
# LIGER can return fewer cells than the concatenated data; 04.12-04.17 keep the overlap
print(sprintf('H cells: %d of %d in the concatenated data (%d not in it)',
              sum(rownames(H_raw) %in% colnames(data)), ncol(data),
              sum(!rownames(H_raw) %in% colnames(data))))
stopifnot(any(rownames(H_raw) %in% colnames(data)))
saveRDS(H, str_glue('{out_path}uinmf_H_raw.rds'))
write.csv(H_raw, str_glue('{out_path}uinmf_H_raw.csv'))

W = getMatrix(data_liger_uinmf$liger, 'W')
saveRDS(W, str_glue('{out_path}uinmf_W.rds'))

# cells every method is scored on (05 / utils.align): in the concatenated data
# with a finite, non-zero H row, i.e. the cells the UINMF methods can keep
usable = rownames(H_raw)[complete.cases(H_raw) & rowSums(abs(H_raw)) > 0]
eval_cells = intersect(colnames(data), usable)
writeLines(eval_cells, str_glue('{out_path}eval_cells.txt'))
print(str_glue('{d}: {length(eval_cells)} / {ncol(data)} cells kept for evaluation'))

liger_uinmf = data_liger_uinmf$seurat
liger_uinmf[['Xenium']] = liger_uinmf[['RNA']]
DefaultAssay(liger_uinmf) = 'Xenium'
liger_uinmf[['RNA']] = NULL

liger_uinmf = liger_uinmf |> JoinLayers()
as.anndata(liger_uinmf, file_path = out_path,
           file_name = '04.11_LIGER_UINMF.h5ad',
           assay = 'Xenium', main_layer = 'counts',
           other_layers = NULL)

print('LIGER UINMF done!')}
