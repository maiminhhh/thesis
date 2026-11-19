setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'bone_marrow'))
source('../library.R')
source('../source.R')
# human
human_dir = 'data/00_raw_data/output-XETG00120_0005071_normal_bone_marrow_20240119_232334'
human = LoadXenium(human_dir,
                    fov = 'fov')
human_cells = read_parquet(file.path(human_dir, 'cells.parquet'))
human$area = human_cells$cell_area[match(colnames(human), human_cells$cell_id)]
stopifnot(!anyNA(human$area))
human_keep_cells = bg_filter(nCount_probe_controls = human$nCount_ControlProbe,
                                nCount_transcripts = human$nCount_Xenium,
                                cell_areas = human$area,
                                n_genes = dim(human)[1],
                                n_total_probes = 20)
human = human[,human_keep_cells]
human_panel = fromJSON(file.path(human_dir, 'gene_panel.json'),
                        simplifyVector = FALSE)
human_targets = human_panel$payload$targets
human_genes = data.frame(
    gene_id = vapply(human_targets, function(t) get_or_na(t$type$data$id), character(1)),
    gene_name = vapply(human_targets, function(t) get_or_na(t$type$data$name), character(1)),
    type = vapply(human_targets, function(t) get_or_na(t$type$descriptor), character(1))
)
human_genes = human_genes |> 
    filter(!is.na(gene_id))
# gene_panel.json is not in the same order as the count matrix: align by name
human_genes = human_genes[match(rownames(human[['Xenium']]), human_genes$gene_name), ]
rownames(human_genes) = NULL
stopifnot(identical(human_genes$gene_name, rownames(human[['Xenium']])))
human[['Xenium']]@meta.data = human_genes
SaveSeuratRds(human, 'data/01_processed_data/01_human_processed.rds')

# mouse
mouse_dir = 'data/00_raw_data/output-XETG00120__0005224__p8_femur__20231122__195355'
mouse = LoadXenium(mouse_dir, 
                    fov = 'fov')
mouse_cells = read_parquet(file.path(mouse_dir, 'cells.parquet'))
mouse$area = mouse_cells$cell_area[match(colnames(mouse), mouse_cells$cell_id)]
stopifnot(!anyNA(mouse$area))
mouse_keep_cells = bg_filter(nCount_probe_controls = mouse$nCount_ControlProbe,
                                nCount_transcripts = mouse$nCount_Xenium,
                                cell_areas = mouse$area,
                                n_genes = dim(mouse)[1],
                                n_total_probes = 20)
mouse = mouse[,mouse_keep_cells]
mouse_panel = fromJSON(file.path(mouse_dir, 'gene_panel.json'),
                        simplifyVector = FALSE)
mouse_targets = mouse_panel$payload$targets
mouse_genes = data.frame(
    gene_id = vapply(mouse_targets, function(t) get_or_na(t$type$data$id), character(1)),
    gene_name = vapply(mouse_targets, function(t) get_or_na(t$type$data$name), character(1)),
    type = vapply(mouse_targets, function(t) get_or_na(t$type$descriptor), character(1))
)
mouse_genes = mouse_genes |> 
    filter(!is.na(gene_id))
# gene_panel.json is not in the same order as the count matrix: align by name
mouse_genes = mouse_genes[match(rownames(mouse[['Xenium']]), mouse_genes$gene_name), ]
rownames(mouse_genes) = NULL
stopifnot(identical(mouse_genes$gene_name, rownames(mouse[['Xenium']])))
mouse[['Xenium']]@meta.data = mouse_genes
SaveSeuratRds(mouse, 'data/01_processed_data/01_mouse_processed.rds')