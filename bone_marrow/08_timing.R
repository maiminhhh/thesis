# Runtime of one R integration method, integration step only (no UMAP, no saving):
#   Rscript 08_timing.R <method> <dataset> <repeat>
# Inputs are read before the clock starts. UINMF and scVI-latent variants read their
# input (UINMF H from 04.11, scVI latent from 04.18) from disk: the time to make that
# input is measured separately (LIGER UINMF here, scVI latent in 08_timing.py).
# Appends one row to data/08_timing/timing.csv (written by run_timing.sh).
source('../library.R')
source('../source.R')
source('integration_sets.R')
options(future.globals.maxSize = 20 * 1024^3)
args = commandArgs(trailingOnly = TRUE)
method = args[1]; s = args[2]; rep = as.integer(args[3])
out_path = integration_sets[[s]]$out

data = readRDS(integration_sets[[s]]$concat)
H = W = NULL
if (grepl(' UINMF$', method)) {
  H = readRDS(str_glue('{out_path}uinmf_H_raw.rds'))
  W = readRDS(str_glue('{out_path}uinmf_W.rds'))
}
if (grepl(' scVI-latent$', method)) {
  H = read_scvi_latent(out_path)
  W = if (grepl('^(CCA|RPCA) ', method))
    matrix(0, nrow(data), ncol(H$hsapiens), dimnames = list(rownames(data), NULL)) else
    matrix(numeric(0), 0, 0)
}
if (method == 'LIGER UINMF') {
  # same preprocessing as 04.11 (outside the clock): shared genes keyed by human Ensembl ID
  human = readRDS(integration_sets[[s]]$human); human$orig.ident = 'hsapiens'
  mouse = readRDS(integration_sets[[s]]$mouse); mouse$orig.ident = 'mmusculus'
  map = data[['Xenium']]@meta.data |>
    dplyr::select(mmusculus_homolog_ensembl_gene, hsapiens_homolog_ensembl_gene)
  mouse[['Xenium']]@meta.data = mouse[['Xenium']]@meta.data |>
    left_join(map, by = join_by(gene_id == mmusculus_homolog_ensembl_gene))
  rownames(mouse) = coalesce(mouse[['Xenium']]@meta.data$hsapiens_homolog_ensembl_gene,
                             mouse[['Xenium']]@meta.data$gene_id) |>
    coalesce(mouse[['Xenium']]@meta.data$gene_name)
  rownames(human) = human[['Xenium']]@meta.data$gene_id |>
    coalesce(human[['Xenium']]@meta.data$gene_name)
}
base = sub(' (UINMF|scVI-latent)$', '', method)
gc()

start = Sys.time()
invisible(switch(if (method == 'LIGER UINMF') method else base,
  'FastMNN'     = integrate_fastmnn(data, H = H, W = W),
  'CCA'         = integrate_cca(data, H = H, W = W),
  'Harmony'     = integrate_harmony(data, H = H, W = W),
  'RPCA'        = integrate_rpca(data, H = H, W = W),
  'LIGER'       = integrate_liger(data),
  'LIGER UINMF' = integrate_liger_uinmf(human, mouse),
  stop('unknown method: ', method)))
sec = as.numeric(difftime(Sys.time(), start, units = 'secs'))

row = data.frame(method = method, dataset = s, rep = rep, seconds = round(sec, 2),
                 device = 'cpu', threads = default_threads(), n_cells = ncol(data),
                 load1 = as.numeric(strsplit(readLines('/proc/loadavg'), ' ')[[1]][1]),
                 time = format(Sys.time(), '%Y-%m-%d %H:%M:%S'))
f = 'data/08_timing/timing.csv'
dir.create(dirname(f), showWarnings = FALSE)
write.table(row, f, sep = ',', row.names = FALSE, col.names = !file.exists(f), append = TRUE)
print(row)
