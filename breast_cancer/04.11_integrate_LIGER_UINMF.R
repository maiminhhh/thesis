# print('1. LIGER UINMF Set up')
# print('--------------------------------------')
# source('../library.R')
# source('../source.R')
# options(future.globals.maxSize = 20 * 1024^3)
# data = readRDS('data/03_concatenated_data.rds')

# human = readRDS('data/02_annotated_data/02_human_annotated.rds')
# human$orig.ident = 'hsapiens'
# mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
# mouse$orig.ident = 'mmusculus'

# out_path = 'data/04_integrated_data/full_mouse/'
# param = MulticoreParam(workers = 64)

# print('2. LIGER UINMF Preprocessing')
# print('--------------------------------------')
# orthologs_mapping = data@assays$MERFISH@meta.data |> 
#   dplyr::select(mmusculus_homolog_ensembl_gene, hsapiens_homolog_ensembl_gene)

# mouse@assays$MERFISH@meta.data = mouse@assays$MERFISH@meta.data |> 
#   left_join(orthologs_mapping, 
#             by = join_by(mmusculus_homolog_ensembl_gene))

# mouse@assays$MERFISH@meta.data$rownames = 
#   coalesce(mouse@assays$MERFISH@meta.data$hsapiens_homolog_ensembl_gene,
#            mouse@assays$MERFISH@meta.data$mmusculus_homolog_ensembl_gene) |> 
#   coalesce(mouse@assays$MERFISH@meta.data$gene_name)

# rownames(mouse) = mouse@assays$MERFISH@meta.data$rownames
# rownames(human) = human[['MERFISH']]@meta.data$ensembl_gene_id |> 
#   coalesce(human[['MERFISH']]@meta.data$gene_name)

# print('3. Run LIGER UINMF')
# print('--------------------------------------')
# data_liger_uinmf = integrate_liger_uinmf(human, mouse)

# print('4. LIGER UINMF Save data')
# print('--------------------------------------')
# saveRDS(data_liger_uinmf$liger, str_glue('{out_path}04.11_LIGER_UINMF.rds'))

# liger_uinmf = data_liger_uinmf$seurat
# liger_uinmf[['MERFISH']] = liger_uinmf[['RNA']]
# DefaultAssay(liger_uinmf) = 'MERFISH'
# liger_uinmf[['RNA']] = NULL

# #write_h5ad(data_liger_uinmf$seurat, str_glue('{out_path}04.11_LIGER_UINMF.h5ad'))
# liger_uinmf = liger_uinmf |> JoinLayers()
# as.anndata(liger_uinmf, file_path = out_path, 
#            file_name = '04.11_LIGER_UINMF.h5ad', 
#            assay = 'MERFISH', main_layer = 'counts', 
#            other_layers = NULL)

# print('5. LIGER UINMF Extract H matrix')
# print('--------------------------------------')
# H = getMatrix(data_liger_uinmf$liger, 'H')
# H = lapply(H, function(M) l2(t(M)))
# H_raw = rbind(H$mmusculus,
#           H$hsapiens)
# saveRDS(H, str_glue('{out_path}uinmf_H_raw.rds'))
# write.csv(H_raw, str_glue('{out_path}uinmf_H_raw.csv'))

# W = getMatrix(data_liger_uinmf$liger, 'W')
# saveRDS(W, str_glue('{out_path}uinmf_W.rds'))

# print('LIGER UINMF done!')

print('1. LIGER UINMF Set up')
print('--------------------------------------')
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
options(future.globals.maxSize = 20 * 1024^3)

# datasets (integration_sets.R): all by default, or e.g. `Rscript 04.11_integrate_LIGER_UINMF.R genes8` to run one dataset (datasets can run in parallel)
source('integration_sets.R')
for (i in run_sets) {
data = readRDS(str_glue('data/{i}/03_concatenated_data.rds'))

human = readRDS(str_glue('data/{i}/02_human_annotated.rds'))
human$orig.ident = 'hsapiens'

mouse = readRDS(str_glue('data/{i}/02_mouse_annotated.rds'))
mouse$orig.ident = 'mmusculus'

out_path = str_glue('data/{i}/')
param = MulticoreParam(workers = 32)

print('2. LIGER UINMF Preprocessing')
print('--------------------------------------')
orthologs_mapping = data[['MERFISH']]@meta.data |> 
  dplyr::select(mmusculus_homolog_ensembl_gene, hsapiens_homolog_ensembl_gene)

mouse[['MERFISH']]@meta.data = mouse[['MERFISH']]@meta.data |> 
  left_join(orthologs_mapping, 
            by = join_by(mmusculus_homolog_ensembl_gene))

mouse[['MERFISH']]@meta.data$rownames = 
  coalesce(mouse@assays$MERFISH@meta.data$hsapiens_homolog_ensembl_gene,
           mouse@assays$MERFISH@meta.data$mmusculus_homolog_ensembl_gene) |> 
  coalesce(mouse@assays$MERFISH@meta.data$gene_name)

rownames(mouse) = mouse@assays$MERFISH@meta.data$rownames
rownames(human) = human[['MERFISH']]@meta.data$ensembl_gene_id |> 
  coalesce(human[['MERFISH']]@meta.data$gene_name)

VariableFeatures(data) = rownames(data)
VariableFeatures(human) = rownames(human)
VariableFeatures(mouse) = rownames(mouse)

print('3. Run LIGER UINMF')
print('--------------------------------------')
data_liger_uinmf = integrate_liger_uinmf(human, mouse)

print('4. LIGER UINMF Save data')
print('--------------------------------------')
saveRDS(data_liger_uinmf$liger, str_glue('{out_path}04.11_LIGER_UINMF.rds'))

liger_uinmf = data_liger_uinmf$seurat
liger_uinmf[['MERFISH']] = liger_uinmf[['RNA']]
DefaultAssay(liger_uinmf) = 'MERFISH'
liger_uinmf[['RNA']] = NULL

#write_h5ad(data_liger_uinmf$seurat, str_glue('{out_path}04.11_LIGER_UINMF.h5ad'))
liger_uinmf = liger_uinmf |> JoinLayers()
as.anndata(liger_uinmf, file_path = out_path, 
           file_name = '04.11_LIGER_UINMF.h5ad', 
           assay = 'MERFISH', main_layer = 'counts', 
           other_layers = NULL)

print('5. LIGER UINMF Extract H matrix')
print('--------------------------------------')
H = getMatrix(data_liger_uinmf$liger, 'H')
H = lapply(H, function(M) l2(t(M)))
# keep H rows with >= 2 non-zero factors (RPCA UINMF gives NA otherwise); every UINMF
# method and the evaluation (eval_cells.txt) must use the same cells
H = lapply(H, function(h) h[rowSums(h != 0) >= 2, , drop = FALSE])
H_raw = rbind(H$mmusculus,
          H$hsapiens)
saveRDS(H, str_glue('{out_path}uinmf_H_raw.rds'))
write.csv(H_raw, str_glue('{out_path}uinmf_H_raw.csv'))

# evaluation cells: cells of the concatenated data with a usable H row
# (LIGER also drops cells with no counts); every 05 script scores on exactly these
eval_cells = intersect(colnames(data), rownames(H_raw))
writeLines(eval_cells, str_glue('{out_path}eval_cells.txt'))
print(str_glue('{i}: {length(eval_cells)} / {ncol(data)} cells kept for evaluation'))

W = getMatrix(data_liger_uinmf$liger, 'W')
saveRDS(W, str_glue('{out_path}uinmf_W.rds'))

print('LIGER UINMF done!')
}