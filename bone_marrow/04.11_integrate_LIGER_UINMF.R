print('1. LIGER UINMF Set up')
print('--------------------------------------')
source('../library.R')
source('../source.R')
source('integration_sets.R')
options(future.globals.maxSize = 20 * 1024^3)
param = MulticoreParam(workers = 32)

for (s in run_sets) {
data = readRDS(integration_sets[[s]]$concat)
human = readRDS(integration_sets[[s]]$human)
human$orig.ident = 'hsapiens'
mouse = readRDS(integration_sets[[s]]$mouse)
mouse$orig.ident = 'mmusculus'
out_path = integration_sets[[s]]$out

print(str_glue('2. LIGER UINMF Preprocessing - {s}'))
print('--------------------------------------')
# shared genes take the human Ensembl ID in both species; unshared keep their own
orthologs_mapping = data[['Xenium']]@meta.data |>
  dplyr::select(mmusculus_homolog_ensembl_gene, hsapiens_homolog_ensembl_gene)

mouse[['Xenium']]@meta.data = mouse[['Xenium']]@meta.data |>
  left_join(orthologs_mapping,
            by = join_by(gene_id == mmusculus_homolog_ensembl_gene))

mouse[['Xenium']]@meta.data$rownames =
  coalesce(mouse[['Xenium']]@meta.data$hsapiens_homolog_ensembl_gene,
           mouse[['Xenium']]@meta.data$gene_id) |>
  coalesce(mouse[['Xenium']]@meta.data$gene_name)

rownames(mouse) = mouse[['Xenium']]@meta.data$rownames
rownames(human) = human[['Xenium']]@meta.data$gene_id |>
  coalesce(human[['Xenium']]@meta.data$gene_name)

print('3. Run LIGER UINMF')
print('--------------------------------------')
data_liger_uinmf = integrate_liger_uinmf(human, mouse)

print('4. LIGER UINMF Save data')
print('--------------------------------------')
saveRDS(data_liger_uinmf$liger, str_glue('{out_path}04.11_LIGER_UINMF.rds'))

liger_uinmf = data_liger_uinmf$seurat
liger_uinmf[['Xenium']] = liger_uinmf[['RNA']]
DefaultAssay(liger_uinmf) = 'Xenium'
liger_uinmf[['RNA']] = NULL
liger_uinmf = liger_uinmf |> JoinLayers()
as.anndata(liger_uinmf, file_path = out_path,
           file_name = '04.11_LIGER_UINMF.h5ad',
           assay = 'Xenium', main_layer = 'counts',
           other_layers = NULL)

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
saveRDS(H, str_glue('{out_path}uinmf_H_raw.rds'))
write.csv(H_raw, str_glue('{out_path}uinmf_H_raw.csv'))

W = getMatrix(data_liger_uinmf$liger, 'W')
saveRDS(W, str_glue('{out_path}uinmf_W.rds'))

# cells every method is scored on (05 / utils.align): in the concatenated data
# with a finite, non-zero H row, i.e. the cells the UINMF methods can keep
usable = rownames(H_raw)[complete.cases(H_raw) & rowSums(abs(H_raw)) > 0]
eval_cells = intersect(colnames(data), usable)
writeLines(eval_cells, str_glue('{out_path}eval_cells.txt'))
print(str_glue('{s}: {length(eval_cells)} / {ncol(data)} cells kept for evaluation'))

print(str_glue('LIGER UINMF done! - {s}'))
}
