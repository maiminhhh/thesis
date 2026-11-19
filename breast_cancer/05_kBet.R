source('../library.R')
source('../source.R')
library(reticulate)
use_condaenv("base", required = TRUE)
scib_m = import("scib_metrics")
nn_module = import("scib_metrics.nearest_neighbors")
all_kbet = list()
## full 
all_kbet$full = data.frame(
    'Algorithm' = character(length = 19),
    'Global kBET' = numeric(length = 19),
    'Per celltype kBET' = numeric(length = 19)
)

path = 'data/04_integrated_data/full_mouse'
f = list.files(path, full.names = TRUE)
f_h5ad = str_detect(f, '.h5ad')
f_rds = str_detect(f, '.rds')
for (i in f[f_h5ad]){
    f_id = str_split(i, '.h5ad')[[1]][1]
    i_rds = str_glue('{f_id}.rds')
    if (!file.exists(i_rds)){
        print(i)
        data = zellkonverter::readH5AD(i, reader = 'R')
        assay = SummarizedExperiment::assayNames(data)
        counts = if ("counts" %in% assay) "counts" else NULL
        norms = intersect(c("data", "logcounts", "X_disp", "X"), assay)[1]
        if (str_detect(str_to_lower(i), 'samap')) {
            rownames(data) = gsub("_", "-", rownames(data))
        }
        data = as.Seurat(data, counts = counts, data = norms)
        SaveSeuratRds(data, i_rds)
    }
}

f = list.files(path, full.names = TRUE)
f = f[str_detect(f, '.rds') & !str_detect(f, 'uinmf_H') & !str_detect(f, 'uinmf_W')]
f = c(f, 'data/03_concatenated_data.rds',
         'data/04_integrated_data/full_mouse/04.14_Harmony_UINMF.rds')
raw = readRDS('data/03_concatenated_data.rds')
out_path = 'data/05_evaluated_score/full_mouse/kbet.rds'
embeddings = c(
    'mnn', NA, 'cca', 'harmony', 'X_scanorama', 'rpca', 'X_scvi', 'X_scanvi', NA, 'inmfNorm', 
    'inmf', 'mnn_uinmf', 'cca_uinmf', 'harmony', 'X_scanorama_uinmf', 'rpca_uinmf', NA,
    'pca', 'uinmf_H'
)
algo = c(
    'FastMNN', 'SAMAP', 'CCA', 'Harmony', 'Scanorama', 'RPCA', 'scVI', 'scANVI', 'BBKNN', 'LIGER',
    'LIGER UINMF', 'FastMNN UINMF', 'CCA UINMF', 'Harmony UINMF', 'Scanorama UINMF', 'RPCA UINMF', 'BBKNN UINMF',
    'Raw', 'Raw UINMF'
)

for (i in 1:length(f)) {
    data = readRDS(f[i])
    name = tools::file_path_sans_ext(basename(f[i]))
    if (str_to_lower(name) |> str_detect('liger_uinmf')) {
        data = data  |> 
            ligerToSeurat(data)
        data@meta.data = raw@meta.data
    }
    if (is.na(embeddings[i])){
        print(f[i])
        all_kbet$full[i,] = c(algo[i], NA, NA)
        next
    }
    if ('celltype' %in% colnames(data@meta.data)){cell_labels = data$celltype}
    else if ('ds_label' %in% colnames(data@meta.data)){cell_labels = data$ds_label}
    else if ('celltype' %in% colnames(raw@meta.data)) {cell_labels = raw@celltype}
    else {cell_labels = raw@ds_label}
    emb = Embeddings(data, reduction = embeddings[i])[,1:20]
    batch_labels = data$orig.ident
    knn_graph = nn_module$pynndescent(emb, n_neighbors = as.integer(100))
    kbet_output = scib_m$kbet(X = knn_graph, batches = batch_labels)
    acceptance_rate = kbet_output[[1]]
    kbet_per_celltype = scib_m$kbet_per_label(X = knn_graph, batches = batch_labels, labels = cell_labels)
    all_kbet$full[i,] = c(algo[i], acceptance_rate, kbet_per_celltype)
}
SaveRDS(all_kbet, 'data/05_evaluated_score/kbet.rds')

# -------------------------------------------------------------------------------------------------------------------------

## dissected 
all_kbet$dissected = data.frame(
    'Algorithm' = character(length = 19),
    'Global kBET' = numeric(length = 19),
    'Per celltype kBET' = numeric(length = 19)
)

path = 'data/04_integrated_data/dissected_mouse'
f = list.files(path, full.names = TRUE)
f_h5ad = str_detect(f, '.h5ad')
f_rds = str_detect(f, '.rds')
for (i in f[f_h5ad]){
    f_id = str_split(i, '.h5ad')[[1]][1]
    i_rds = str_glue('{f_id}.rds')
    if (!file.exists(i_rds)){
        print(i)
        data = zellkonverter::readH5AD(i, reader = 'R')
        assay = SummarizedExperiment::assayNames(data)
        counts = if ("counts" %in% assay) "counts" else NULL
        norms = intersect(c("data", "logcounts", "X_disp", "X"), assay)[1]
        if (str_detect(str_to_lower(i), 'samap')) {
            rownames(data) = gsub("_", "-", rownames(data))
        }
        data = as.Seurat(data, counts = counts, data = norms)
        SaveSeuratRds(data, i_rds)
    }
}

f = list.files(path, full.names = TRUE)
f = f[str_detect(f, '.rds') & !str_detect(f, 'uinmf_H') & !str_detect(f, 'uinmf_W')]
raw = readRDS('data/03_concatenated_data_subset.rds')
f = c(f, 'data/03_concatenated_data_subset.rds',
         'data/04_integrated_data/dissected_mouse/04.14_Harmony_UINMF.rds')
out_path = 'data/05_evaluated_score/dissected_mouse/kbet.rds'
embeddings = c(
    'mnn', NA, 'cca', 'harmony', 'X_scanorama', 'rpca', 'X_scvi', 'X_scanvi', NA, 'inmfNorm', 
    'inmf', 'mnn_uinmf', 'cca_uinmf', 'harmony', 'X_scanorama_uinmf', 'rpca_uinmf', NA,
    'pca', 'uinmf_H'
)
algo = c(
    'FastMNN', 'SAMAP', 'CCA', 'Harmony', 'Scanorama', 'RPCA', 'scVI', 'scANVI', 'BBKNN', 'LIGER',
    'LIGER UINMF', 'FastMNN UINMF', 'CCA UINMF', 'Harmony UINMF', 'Scanorama UINMF', 'RPCA UINMF', 'BBKNN UINMF',
    'Raw', 'Raw UINMF'
)

for (i in 1:length(f)) {
    data = readRDS(f[i])
    name = tools::file_path_sans_ext(basename(f[i]))
    if (str_to_lower(name) |> str_detect('liger_uinmf')) {
        data = data  |> 
            ligerToSeurat(data)
        data@meta.data = raw@meta.data
    }
    if (is.na(embeddings[i])){
        print(f[i])
        all_kbet$dissected[i,] = c(algo[i], NA, NA)
        next
    }
    if ('celltype' %in% colnames(data@meta.data)){cell_labels = data$celltype}
    else if ('ds_label' %in% colnames(data@meta.data)){cell_labels = data$ds_label}
    else if ('celltype' %in% colnames(raw@meta.data)) {cell_labels = raw@celltype}
    else {cell_labels = raw@ds_label}
    emb = Embeddings(data, reduction = embeddings[i])[,1:20]
    batch_labels = data$orig.ident
    knn_graph = nn_module$pynndescent(emb, n_neighbors = as.integer(100))
    kbet_output = scib_m$kbet(X = knn_graph, batches = batch_labels)
    acceptance_rate = kbet_output[[1]]
    kbet_per_celltype = scib_m$kbet_per_label(X = knn_graph, batches = batch_labels, labels = cell_labels)
    all_kbet$dissected[i,] = c(algo[i], acceptance_rate, kbet_per_celltype)
}
SaveRDS(all_kbet, 'data/05_evaluated_score/kbet.rds')