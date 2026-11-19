# Rerun the annotators that use the human (CELLxGENE) reference, after the expm1(X) fix:
#   human: rebuild human_processed_human_use.rds (02_human_annotate.qmd), then celltype_seurat / celltype_singleR
#   mouse: rebuild human_processed_mouse_use.rds, then celltype_seurat_cross / celltype_singleR_cross
# The code is taken from the 02 files by content, so it is exactly what the scripts do.
# usage: Rscript reannotate.R human|mouse
which = commandArgs(TRUE)[1]
setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))
pdf(NULL)

# lines of `file` from the line equal to `from` to the line equal to `to` (both trimmed, each unique)
code_between = function(file, from, to) {
  src = readLines(file)
  i = which(trimws(src) == from); j = which(trimws(src) == to)
  if (length(i) != 1 || length(j) != 1 || j < i) stop(file, ': cannot find a unique "', from, '" ... "', to, '"')
  src[i:j]
}
run = function(code, label) {
  message(sprintf('\n===== %s [%s] =====', label, format(Sys.time())))
  eval(parse(text = code), envir = globalenv())
}

run(code_between('02_human_annotated.R', "setwd(file.path(Sys.getenv('PROJECT_ROOT', '..'), 'colon'))", "source('../source.R')"), 'setup')

if (which == 'human') {
  human = readRDS('data/02_annotated_data/02_human_annotated.rds')
  run(code_between('02_human_annotate.qmd', "human_ref = read_h5ad('data/02_annotation_refs/human.h5ad')",
                   "SaveSeuratRds(human_ref, 'data/02_annotation_refs/human_processed_human_use.rds')"),
      'build human reference (02_human_annotate.qmd)')
  human = readRDS('data/02_annotated_data/02_human_annotated.rds')
  rownames(human) = human@assays$Xenium@meta.data$ensembl_id
  old = human@meta.data[c('celltype_seurat', 'celltype_singleR')]
  run(code_between('02_human_annotated.R', "human_ref = readRDS('data/02_annotation_refs/human_processed_human_use.rds')",
                   "print('Done SingleR!')"),
      'same-species Seurat + SingleR (02_human_annotated.R)')
  new = human@meta.data[c('celltype_seurat', 'celltype_singleR')]
} else {
  run(code_between('02_mouse_annotate.qmd', "mouse_ref = read_h5ad('data/02_annotation_refs/human.h5ad')",
                   "SaveSeuratRds(mouse_ref_mouse, 'data/02_annotation_refs/human_processed_mouse_use.rds')"),
      'build cross-species reference (02_mouse_annotate.qmd)')
  mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
  old = mouse@meta.data[c('celltype_seurat_cross', 'celltype_singleR_cross')]
  src = readLines('02_mouse_annotate.qmd')
  # the "## Seurats" and "## SingleR" chunks of the cross-species section
  seurat_start = which(trimws(src) == "mouse_seurat_label = LabelTransfer(data = mouse, ref = mouse_ref)")
  singler_start = which(trimws(src) == "mouse_singleR_label = LabelSingleR(data = mouse, ref = mouse_ref)")
  stopifnot(length(seurat_start) == 2, length(singler_start) == 2)   # within-species, then cross
  ends = which(trimws(src) == "SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')")
  mouse_ref = readRDS('data/02_annotation_refs/human_processed_mouse_use.rds')
  rownames(mouse) = mouse@assays$Xenium@meta.data$ensembl_id
  run(src[seurat_start[2]:min(ends[ends > seurat_start[2]])], 'cross-species Seurat (02_mouse_annotate.qmd)')
  run(src[singler_start[2]:min(ends[ends > singler_start[2]])], 'cross-species SingleR (02_mouse_annotate.qmd)')
  stopifnot(grepl('cross', paste(src[seurat_start[2]:min(ends[ends > seurat_start[2]])], collapse = ' ')),
            grepl('cross', paste(src[singler_start[2]:min(ends[ends > singler_start[2]])], collapse = ' ')))
  mouse = readRDS('data/02_annotated_data/02_mouse_annotated.rds')
  new = mouse@meta.data[c('celltype_seurat_cross', 'celltype_singleR_cross')]
}

# how much the fix changed each annotator
for (a in names(old)) {
  message(sprintf('\n%s: %.1f%% of cells changed label', a,
                  100 * mean(as.character(old[[a]]) != as.character(new[[a]]), na.rm = TRUE)))
  print(table(old = old[[a]], new = new[[a]], useNA = 'ifany'))
}
message('reannotation done [', format(Sys.time()), ']')
