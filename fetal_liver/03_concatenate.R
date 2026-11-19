setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'fetal_liver'))
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))   # anndata, for as.anndata()
source('../library.R')
source('../source.R')
options(Seurat.nthreads = 16)   # uwot UMAP (Seurat defaults to 1 thread)

human = readRDS('data/02_annotated_data/02_human_annotated.rds')
mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
homology_tbl = read.csv('../homology_tbl.csv')

# dissect the mouse section by its coordinates (window chosen in visualisations.qmd)
xmin = 0;  xmax = 3100
ymin = 500; ymax = 3000
coords = GetTissueCoordinates(mouse)
keep_cells = coords$cell[coords$x >= xmin & coords$x <= xmax &
                         coords$y >= ymin & coords$y <= ymax]
mouse_subset = mouse |>
  subset(cells = keep_cells)
SaveSeuratRds(mouse_subset, 'data/02_annotated_data/02_mouse_subset_annotated.rds')
as.anndata(mouse_subset, file_path = 'data/02_annotated_data/',
           file_name = '02_mouse_subset_annotated.h5ad', assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))
print(sprintf('Dissected mouse: %d of %d cells', ncol(mouse_subset), ncol(mouse)))

rownames(mouse) = mouse[['Xenium']]@meta.data$gene_id
rownames(mouse_subset) = mouse_subset[['Xenium']]@meta.data$gene_id

rownames(human) = human[['Xenium']]@meta.data$gene_id

# per-species inputs (Ensembl ids, for SAMap / LIGER UINMF) restricted to the
# cells kept in the concatenated data, saved as 02_{human,mouse}_annotated in out_dir
save_species_inputs = function(concat, human, mouse, out_dir) {
  # concat cell names are '<species>_<cell id>' (add.cell.ids in concatenate)
  kept = colnames(concat)
  human_cells = sub('^hsapiens_', '', kept[concat$orig.ident == 'hsapiens'])
  mouse_cells = sub('^mmusculus_', '', kept[concat$orig.ident == 'mmusculus'])
  stopifnot(all(human_cells %in% colnames(human)), all(mouse_cells %in% colnames(mouse)))
  human = human[, human_cells]
  mouse = mouse[, mouse_cells]
  stopifnot(ncol(human) + ncol(mouse) == ncol(concat))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  SaveSeuratRds(human, file.path(out_dir, '02_human_annotated.rds'))
  SaveSeuratRds(mouse, file.path(out_dir, '02_mouse_annotated.rds'))
  as.anndata(human, file_path = out_dir,
             file_name = '02_human_annotated.h5ad', assay = 'Xenium',
             main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')
  as.anndata(mouse, file_path = out_dir,
             file_name = '02_mouse_annotated.h5ad', assay = 'Xenium',
             main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')
}

# full section
concat = concatenate(human, mouse, homology_tbl)
VariableFeatures(concat) = rownames(concat)
concat = JoinLayers(concat)
concat$nCount_Xenium = colSums(x = concat[['Xenium']]$counts)
concat = concat |>
  subset(nCount_Xenium > 0)

concat = concat |>
  RunPCA() |>
  RunUMAP(dims = 1:20)

SaveSeuratRds(concat, 'data/03_concatenated_data/03_concatenated_data_full.rds')
as.anndata(concat, file_path = 'data/03_concatenated_data/',
           file_name = '03_concatenated_data_full.h5ad', assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))
save_species_inputs(concat, human, mouse, 'data/03_concatenated_data/full_mouse/')
# subsection
concat = concatenate(human, mouse_subset, homology_tbl)
VariableFeatures(concat) = rownames(concat)
concat = JoinLayers(concat)
concat$nCount_Xenium = colSums(x = concat[['Xenium']]$counts)
concat = concat |>
  subset(nCount_Xenium > 0)

concat = concat |>
  RunPCA() |>
  RunUMAP(dims = 1:20)

SaveSeuratRds(concat, 'data/03_concatenated_data/03_concatenated_data_subset.rds')
as.anndata(concat, file_path = 'data/03_concatenated_data/',
           file_name = '03_concatenated_data_subset.h5ad', assay = 'Xenium',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))
save_species_inputs(concat, human, mouse_subset, 'data/03_concatenated_data/dissected_mouse/')

# subset genes
# p = fraction of shared genes DROPPED; saved to data/genes<10p>/
mouse = mouse_subset
shared_genes = concat[['Xenium']][[]]
# mouse id for every shared (human id) gene, looked up by name rather than by row order
hs_to_mm = setNames(shared_genes$mmusculus_homolog_ensembl_gene,
                    shared_genes$hsapiens_homolog_ensembl_gene)
stopifnot(all(rownames(concat) %in% names(hs_to_mm)))
n_shared = nrow(concat)
p = c(0.8, 0.6, 0.4, 0.2)
# nested drops: one random gene order, each level drops a longer prefix of it
set.seed(1)
gene_order = sample(rownames(concat))
# levels run one after another: forking them (mclapply) segfaulted in RunPCA and
# RunUMAP, because OpenMP/OpenBLAS are not safe in a forked child
done = lapply(p, FUN = function(frac) {
    drop = gene_order[seq_len(round(frac * n_shared))]
    # renormalise each species on its remaining genes (log-normalisation and scaling
    # depend on which genes are present), then concatenate them as for full / dissected
    renormalise = function(obj) {
      VariableFeatures(obj) = rownames(obj)
      obj |> NormalizeData(verbose = FALSE) |> ScaleData(verbose = FALSE)
    }
    human_p = renormalise(human[!rownames(human) %in% drop, ])
    mouse_p = renormalise(mouse[!rownames(mouse) %in% hs_to_mm[drop], ])
    stopifnot(nrow(human) - nrow(human_p) == length(drop),
              nrow(mouse) - nrow(mouse_p) == length(drop))
    concat_p = concatenate(human_p, mouse_p, homology_tbl)
    stopifnot(setequal(rownames(concat_p), setdiff(rownames(concat), drop)))
    VariableFeatures(concat_p) = rownames(concat_p)
    concat_p = JoinLayers(concat_p)
    # cells left with no counts on the kept shared genes are removed everywhere
    concat_p$nCount_Xenium = colSums(x = concat_p[['Xenium']]$counts)
    concat_p = concat_p |>
      subset(nCount_Xenium > 0)
    concat_p = concat_p |>
      RunPCA(npcs = min(50, nrow(concat_p) - 1)) |>
      # PCA init: with few genes many cells are identical, the kNN graph falls apart and
      # uwot's spectral init hangs or segfaults in RSpectra (UMAP here is only for viewing)
      RunUMAP(dims = 1:min(20, nrow(concat_p) - 1), uwot.init = 'pca')

    out_dir = sprintf('data/genes%d/', round(10 * frac))
    save_species_inputs(concat_p, human_p, mouse_p, out_dir)
    SaveSeuratRds(concat_p, file.path(out_dir, '03_concatenated_data.rds'))
    as.anndata(concat_p, file_path = out_dir,
               file_name = '03_concatenated_data.h5ad', assay = 'Xenium',
               main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')
    sprintf('Dropped %d of %d shared genes -> %s (cells kept: %d of %d)',
            length(drop), n_shared, out_dir, ncol(concat_p), ncol(concat))
})
for (msg in done) print(msg)
