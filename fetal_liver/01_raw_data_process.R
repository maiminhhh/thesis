setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'fetal_liver'))
source('../library.R')
source('../source.R')

process_xenium = function(data_dir, out_file) {
    obj = LoadXenium(data_dir, fov = 'fov')
    cells = read_parquet(file.path(data_dir, 'cells.parquet'))
    obj$area = cells$cell_area[match(colnames(obj), cells$cell_id)]
    stopifnot(!anyNA(obj$area))
    keep_cells = bg_filter(nCount_probe_controls = obj$nCount_ControlProbe,
                           nCount_transcripts = obj$nCount_Xenium,
                           cell_areas = obj$area,
                           n_genes = nrow(obj[['Xenium']]),
                           n_total_probes = nrow(obj[['ControlProbe']]))
    obj = obj[, keep_cells]
    panel = fromJSON(file.path(data_dir, 'gene_panel.json'),
                     simplifyVector = FALSE)
    targets = panel$payload$targets
    genes = data.frame(
        gene_id = vapply(targets, function(t) get_or_na(t$type$data$id), character(1)),
        gene_name = vapply(targets, function(t) get_or_na(t$type$data$name), character(1)),
        type = vapply(targets, function(t) get_or_na(t$type$descriptor), character(1))
    )
    genes = genes |>
        filter(!is.na(gene_id))
    # gene_panel.json is not in the same order as the count matrix: align by name
    genes = genes[match(rownames(obj[['Xenium']]), genes$gene_name), ]
    rownames(genes) = NULL
    stopifnot(identical(genes$gene_name, rownames(obj[['Xenium']])))
    obj[['Xenium']]@meta.data = genes
    SaveSeuratRds(obj, out_file)
    obj
}

# human
human = process_xenium('data/00_raw_data/output-XETG00120_0005059_fetal_liver_20240119_232334',
                       'data/01_processed_data/01_human_processed.rds')

# mouse
mouse = process_xenium('data/00_raw_data/output-XETG00120__0005224__E145_liver__20231122__195355',
                       'data/01_processed_data/01_mouse_processed.rds')
