source('../library.R')
source('../source.R')
library(reticulate)
use_condaenv("base", required = TRUE)
scib_m = import("scib_metrics")
nn_module = import("scib_metrics.nearest_neighbors")

setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
source('integration_datasets.R')

# one row per method, looked up by file name (not by position in the folder listing);
# from_h5ad = written by Python, so its .rds is converted from the .h5ad;
# embedding = NA for graph-only methods (no embedding to run kBET on)
methods = tribble(
  ~file,                      ~algo,             ~embedding,           ~from_h5ad,
  '04.01_FastMNN',            'FastMNN',         'mnn',                FALSE,
  '04.02_SAMap',              'SAMAP',           NA,                   TRUE,
  '04.03_CCA',                'CCA',             'cca',                FALSE,
  '04.04_Harmony',            'Harmony',         'harmony',            FALSE,
  '04.05_Scanorama',          'Scanorama',       'X_scanorama',        TRUE,
  '04.06_RPCA',               'RPCA',            'rpca',               FALSE,
  '04.07_scVI',               'scVI',            'X_scvi',             TRUE,
  '04.08_scANVI',             'scANVI',          'X_scanvi',           TRUE,
  '04.09_BBKNN',              'BBKNN',           NA,                   TRUE,
  '04.10_LIGER',              'LIGER',           'inmfNorm',           FALSE,
  '04.11_LIGER_UINMF',        'LIGER UINMF',     'inmf',               FALSE,
  '04.12_FastMNN_UINMF',      'FastMNN UINMF',   'mnn_uinmf',          FALSE,
  '04.13_CCA_UINMF',          'CCA UINMF',       'cca_uinmf',          FALSE,
  '04.14_Harmony_UINMF',      'Harmony UINMF',   'harmony',            FALSE,
  '04.15_Scanorama_UINMF',    'Scanorama UINMF', 'X_scanorama_uinmf',  TRUE,
  '04.16_RPCA_UINMF',         'RPCA UINMF',      'rpca_uinmf',         FALSE,
  '04.17_BBKNN_UINMF',        'BBKNN UINMF',     NA,                   TRUE,
  'raw',                      'Raw',             'pca',                FALSE,
  '04.14_Harmony_UINMF',      'Raw UINMF',       'uinmf_H',            FALSE
)

all_kbet = list()
for (d in names(datasets)) {
  path = sub('/$', '', datasets[[d]]$out)
  raw_file = str_glue('{datasets[[d]]$concat}.rds')
  raw = readRDS(raw_file)
  all_kbet[[d]] = data.frame('Algorithm' = methods$algo,
                             'Global kBET' = NA_real_,
                             'Per celltype kBET' = NA_real_,
                             'n cells' = NA_integer_,
                             check.names = FALSE)

  for (i in seq_len(nrow(methods))) {
    m = methods[i, ]
    if (is.na(m$embedding)) {
      print(str_glue('{d}: {m$algo} is graph-only, skipped'))
      next
    }
    f_rds = if (m$file == 'raw') raw_file else str_glue('{path}/{m$file}.rds')

    # Python outputs: (re)convert when the .rds is missing or older than the .h5ad,
    # so a rerun of 04 is never scored from a stale conversion
    if (m$from_h5ad) {
      f_h5ad = str_glue('{path}/{m$file}.h5ad')
      if (!file.exists(f_rds) || file.mtime(f_rds) < file.mtime(f_h5ad)) {
        print(str_glue('converting {f_h5ad}'))
        data = zellkonverter::readH5AD(f_h5ad, reader = 'R')
        assay = SummarizedExperiment::assayNames(data)
        counts = if ("counts" %in% assay) "counts" else NULL
        norms = intersect(c("data", "logcounts", "X_disp", "X"), assay)[1]
        data = as.Seurat(data, counts = counts, data = norms)
        SaveSeuratRds(data, f_rds)
      }
    }

    data = readRDS(f_rds)
    # 04.11 saves the liger object itself
    if (inherits(data, 'liger')) data = ligerToSeurat(data)

    # UINMF-based methods can hold fewer cells: score the cells shared with raw,
    # taking batch and cell type from raw by cell name
    cells = intersect(colnames(data), colnames(raw))
    emb = Embeddings(data, reduction = m$embedding)[cells, , drop = FALSE]
    emb = emb[, seq_len(min(20, ncol(emb))), drop = FALSE]
    batch_labels = as.character(raw@meta.data[cells, 'orig.ident'])
    cell_labels = as.character(raw@meta.data[cells, 'celltype'])
    print(str_glue('{d}: {m$algo} | {length(cells)} / {ncol(raw)} cells | {ncol(emb)} dims'))

    knn_graph = nn_module$pynndescent(emb, n_neighbors = as.integer(100))
    kbet_output = scib_m$kbet(X = knn_graph, batches = batch_labels)
    kbet_per_celltype = scib_m$kbet_per_label(X = knn_graph, batches = batch_labels, labels = cell_labels)
    all_kbet[[d]][i, 'Global kBET'] = kbet_output[[1]]
    all_kbet[[d]][i, 'Per celltype kBET'] = kbet_per_celltype
    all_kbet[[d]][i, 'n cells'] = length(cells)
  }
  print(all_kbet[[d]])
  saveRDS(all_kbet[[d]], str_glue('{datasets[[d]]$scores}05_kbet.rds'))
}
