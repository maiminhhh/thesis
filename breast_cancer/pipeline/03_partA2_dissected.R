# generated from 03_concatenate.qmd: the dissected dataset only, starting from the saved full_mouse mouse object
# (= chunk 1's prefixed, Ensembl-named mouse restricted to the concatenated cells)
Sys.setenv(RETICULATE_PYTHON = file.path(Sys.getenv('CONDA_ENVS', path.expand('~/miniconda3/envs')), 'torch/bin/python'))
source('../library.R')
source('../source.R')
drop_reductions = function(obj) { for (r in Reductions(obj)) obj[[r]] = NULL; obj }
mouse = readRDS('data/full_mouse/02_mouse_annotated.rds') |> drop_reductions()
coords = GetTissueCoordinates(mouse)   # has a `cell` column; rownames are just 1..n
coords$celltype = mouse$celltype[coords$cell]
# dissected mouse = central tumour block (../EBCC_PIPELINE.md section 5b): the human sample is a breast tumour
# TMA core, so the matching mouse compartment is the tumour mass; this rectangle stays inside it and leaves out
# the normal mammary gland, fat and lymph-node-like tissue at low x and the capsule rim at the edges
# (map: ../figures/breast_cancer/spatial/mouse_dissection_map_*.png). Previous cut: 8560-11745 x 920-5455.
# GetTissueCoordinates units; in ImageDimPlot y is the horizontal and x the vertical axis.
xmin = 4500; xmax = 8500
ymin = 1500; ymax = 5500

coords$keep = coords$x >= xmin &
  coords$x <= xmax & 
  coords$y >= ymin &
  coords$y <= ymax

mouse_subset = mouse |> 
  subset(cells = coords$cell[coords$keep])

SaveSeuratRds(mouse_subset, 'data/03_mouse_subset.rds')
as.anndata(mouse_subset, file_path = 'data/',
           file_name = '03_mouse_subset.h5ad', assay = 'MERFISH',
           main_layer = 'counts', other_layers = c('data', 'scale.data'))


source('../library.R')
source('../source.R')
homology_tbl = read.csv('../homology_tbl.csv')
# the 02 PCA/UMAP of the species objects are not used here (03 makes its own) and break gene / cell renames
drop_reductions = function(obj) { for (r in Reductions(obj)) obj[[r]] = NULL; obj }
human = readRDS('data/02_annotated_data/02_human_annotated.rds') |> drop_reductions()
ids = coalesce(human[['MERFISH']]@meta.data$ensembl_gene_id,
               human[['MERFISH']]@meta.data$gene_name)

rownames(human) = ids
mouse_subset = readRDS('data/03_mouse_subset.rds') |> drop_reductions()
# mouse_subset was cut from the already-prefixed mouse; drop it so concatenate() adds it only once
mouse_subset = RenameCells(mouse_subset, new.names = sub('^mmusculus_', '', colnames(mouse_subset)))
rownames(mouse_subset) = coalesce(mouse_subset[['MERFISH']]@meta.data$mmusculus_homolog_ensembl_gene, 
                                  mouse_subset[['MERFISH']]@meta.data$gene_name)
concat = concatenate(human, mouse_subset, homology_tbl,
                     tech_human = 'MERFISH', tech_mouse = 'MERFISH')
VariableFeatures(concat) = rownames(concat)
concat = JoinLayers(concat)
concat$nCount_MERFISH = colSums(x = concat[['MERFISH']]$counts)
concat = concat |> 
  subset(nCount_MERFISH > 0) |> 
  RunPCA(npcs = min(50, nrow(concat) - 1)) |>     # PCA + UMAP of the unintegrated data, used by the figures
  RunUMAP(dims = 1:20)

mouse_subset = RenameCells(mouse_subset, new.names = paste0("mmusculus_", colnames(mouse_subset)))

human = RenameCells(human, new.names = paste0("hsapiens_", colnames(human)))

mouse_keep_cells = intersect(colnames(concat), 
                             colnames(mouse_subset))
human_keep_cells = intersect(colnames(concat),
                             colnames(human))

mouse_subset = mouse_subset[,mouse_keep_cells]
human = human[,human_keep_cells]

# dissected mouse dataset
SaveSeuratRds(mouse_subset, 'data/dissected_mouse/02_mouse_annotated.rds')
SaveSeuratRds(human, 'data/dissected_mouse/02_human_annotated.rds')

as.anndata(mouse_subset, file_path = 'data/dissected_mouse/', 
           file_name = '02_mouse_annotated.h5ad', assay = 'MERFISH', 
           main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')
as.anndata(human, file_path = 'data/dissected_mouse/', 
           file_name = '02_human_annotated.h5ad', assay = 'MERFISH', 
           main_layer = 'counts', other_layers = c('data', 'scale.data'), mode = 'w')

SaveSeuratRds(concat, 'data/dissected_mouse/03_concatenated_data.rds')
as.anndata(concat, file_path = 'data/dissected_mouse/', file_name = '03_concatenated_data.h5ad', 
           assay = 'MERFISH', main_layer = 'counts', 
           other_layers = c('data', 'scale.data'))


