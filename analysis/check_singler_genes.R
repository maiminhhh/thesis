# Does uneven panel coverage go with SingleR's over-assignment? Per dataset (same-species reference), SingleR is
# retrained exactly as in the annotation (source.R LabelSingleR: reference 'data' layer, genes shared with the query,
# de.method = 't', default de.n) to recover the genes it uses to recognise each reference type: the union, over all
# pairwise comparisons, of the genes it picks as up-regulated in that type. Only the reference is used for training;
# the query is not relabelled. Check: classifySingleR with the retrained model on a random 3,000 query cells is
# compared with the saved SingleR labels. Detection: mean over a type's genes of the share of query cells detecting
# the gene. Over-assignment of a type: share of cells SingleR gives it / share Seurat
# gives it (all cells both label; SingleR pruned calls left out).
# Run from the repository root with ~/bin/Rscript45 (R 4.5).
# Usage: Rscript45 analysis/check_singler_genes.R <tissue> <species> -> figures/cross_tissue/singler_genes/<tissue>_<species>.csv
setwd(Sys.getenv('PROJECT_ROOT', '.')); source('library.R')
suppressPackageStartupMessages({ library(Seurat); library(SingleR); library(dplyr); library(readr) })

cfg = list(
  colon_human = list(ref = 'colon/data/02_annotation_refs/human_processed_human_use.rds', col = 'cell_type_regroup', tech = 'Xenium', gene_names = FALSE, id = 'ensembl_id'),
  colon_mouse = list(ref = 'colon/data/02_annotation_refs/mouse_processed_mouse_use.rds', col = 'cell_type_regroup', tech = 'Xenium', gene_names = FALSE, id = 'ensembl_id'),
  bone_marrow_human = list(ref = 'bone_marrow/data/02_annotation_refs/02_human_ref.rds', col = 'celltype', tech = 'Xenium', gene_names = TRUE),
  bone_marrow_mouse = list(ref = 'bone_marrow/data/02_annotation_refs/02_mouse_ref.rds', col = 'celltype', tech = 'Xenium', gene_names = TRUE),
  fetal_liver_human = list(ref = 'fetal_liver/data/02_annotation_refs/02_human_ref.rds', col = 'celltype', tech = 'Xenium', gene_names = TRUE),
  fetal_liver_mouse = list(ref = 'fetal_liver/data/02_annotation_refs/02_mouse_ref.rds', col = 'celltype', tech = 'Xenium', gene_names = TRUE),
  breast_cancer_human = list(ref = 'breast_cancer/data/02_annotation_refs/02_human_ref.rds', col = 'celltype_major', tech = 'MERFISH', gene_names = FALSE))

args = commandArgs(TRUE)
t = args[1]; s = args[2]; k = cfg[[paste(t, s, sep = '_')]]
out_dir = 'figures/cross_tissue/singler_genes'; dir.create(out_dir, showWarnings = FALSE)
obj = readRDS(sprintf('%s/data/02_annotated_data/02_%s_annotated.rds', t, s))
# query genes as SingleR saw them: bone marrow and fetal liver renamed rows to gene names before annotating
if (k$gene_names) rownames(obj) = obj[[k$tech]]@meta.data$gene_name
# colon was annotated while rows were Ensembl IDs (the colon references use them); symbols kept for reporting
symbol = setNames(obj[[k$tech]]@meta.data$gene_name, rownames(obj))
if (!is.null(k$id)) { symbol = setNames(obj[[k$tech]]@meta.data$gene_name, obj[[k$tech]]@meta.data[[k$id]]); rownames(obj) = names(symbol) }
if (inherits(obj[[k$tech]], 'Assay5')) obj[[k$tech]] = JoinLayers(obj[[k$tech]])
Q = LayerData(obj, assay = k$tech, layer = 'data')
meta = obj@meta.data; rm(obj); gc()
ref = readRDS(k$ref)
if (inherits(ref[['RNA']], 'Assay5')) ref[['RNA']] = JoinLayers(ref[['RNA']])
R = LayerData(ref, assay = 'RNA', layer = 'data'); labels = as.character(ref@meta.data[[k$col]]); rm(ref); gc()
common = intersect(rownames(Q), rownames(R))
message(sprintf('%s %s: %d query genes, %d reference genes, %d shared', t, s, nrow(Q), nrow(R), length(common)))
trained = trainSingleR(R[common, ], labels = labels, de.method = 't')
up = lapply(trained$markers$full, \(per_other) unique(unlist(per_other)))
# reproduce check on a random sample of query cells
set.seed(1); idx = sample(ncol(Q), min(3000, ncol(Q)))
pred = classifySingleR(Q[common, idx], trained)
saved = as.character(meta$celltype_singleR)[idx]
ok = !is.na(saved) & !is.na(pred$pruned.labels)
match = mean(pred$pruned.labels[ok] == saved[ok])
message(sprintf('reproduction: %.1f%% of %d sampled cells get the saved SingleR label', 100 * match, sum(ok)))
# shares of cells per label, both annotators, cells both label
both = !is.na(meta$celltype_singleR) & !is.na(meta$celltype_seurat)
sh = \(x) prop.table(table(factor(as.character(x)[both], levels = names(up))))
# how often each type's SingleR genes are detected in the query: per gene, share of query cells with a count > 0;
# per type, the mean over its genes
det = Matrix::rowMeans(Q[common, ] > 0)
res = tibble(tissue = t, species = s, type = names(up), n_singler_genes = lengths(up), n_shared_genes = length(common),
             mean_detect = vapply(up, \(g) mean(det[g]), numeric(1)),
             n_ref_cells = as.integer(table(factor(labels, levels = names(up)))),
             share_singler = as.numeric(sh(meta$celltype_singleR)), share_seurat = as.numeric(sh(meta$celltype_seurat)),
             reproduction = match, genes = vapply(up, \(g) paste(symbol[g], collapse = ';'), character(1)))
write_csv(res, file.path(out_dir, sprintf('%s_%s.csv', t, s)))
print(as.data.frame(select(res, type, n_singler_genes, mean_detect, n_ref_cells, share_singler, share_seurat)), digits = 3)
