# Panel genes detected per cell, on the post-QC data that was annotated (02 annotated objects, counts layer).
# Per tissue and species: panel genes in the object, cells, median (and quartiles) of genes detected per cell
# (count > 0), the same as a share of the panel, and panel genes with zero counts in every post-QC cell.
# Run from the repository root. Reads saved objects only.
# Usage: Rscript analysis/check_genes_per_cell.R <tissue> <species>  -> figures/cross_tissue/genes_per_cell/<tissue>_<species>.csv
#        Rscript analysis/check_genes_per_cell.R                     -> figures/cross_tissue/genes_per_cell.csv (combined), printed
setwd(Sys.getenv('PROJECT_ROOT', '.')); source('library.R')
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(purrr); library(readr) })
setwd(Sys.getenv('PROJECT_ROOT', '.'))
out_dir = 'figures/cross_tissue/genes_per_cell'
args = commandArgs(TRUE)
if (length(args) == 2) {
  dir.create(out_dir, showWarnings = FALSE)
  obj = readRDS(sprintf('%s/data/02_annotated_data/02_%s_annotated.rds', args[1], args[2]))
  X = LayerData(obj, assay = DefaultAssay(obj), layer = 'counts')
  genes = if ('gene_name' %in% names(obj[[DefaultAssay(obj)]]@meta.data)) obj[[DefaultAssay(obj)]]@meta.data$gene_name else rownames(X)
  rm(obj); gc()
  per_cell = Matrix::colSums(X > 0)
  per_gene = Matrix::rowSums(X)
  q = quantile(per_cell, c(0.25, 0.5, 0.75))
  write_csv(tibble(tissue = args[1], species = args[2], n_cells = ncol(X), n_genes = nrow(X),
                   q25_genes = q[[1]], median_genes = q[[2]], q75_genes = q[[3]],
                   median_pct_panel = 100 * q[[2]] / nrow(X), mean_genes = mean(per_cell),
                   zero_genes = sum(per_gene == 0), zero_gene_names = paste(genes[per_gene == 0], collapse = ';')),
            file.path(out_dir, sprintf('%s_%s.csv', args[1], args[2])))
  quit(save = 'no')
}
d = map_dfr(list.files(out_dir, '\\.csv$', full.names = TRUE), read_csv, show_col_types = FALSE)
write_csv(d, 'figures/cross_tissue/genes_per_cell.csv')
print(as.data.frame(d), digits = 3)
