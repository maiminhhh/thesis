# datasets every 04 script runs on (mirrors integration_datasets.py): comment out an entry to skip it
#   concat = concatenated input from 03 (no extension: .rds for R, .h5ad for python)
#   inputs = folder with the per-species 02_{human,mouse}_annotated objects (Ensembl ids),
#            restricted to the cells kept in concat
#   out    = output folder
#   scores = folder the 05 evaluation scores are written to
datasets = list(
  full_mouse      = list(concat = 'data/03_concatenated_data/03_concatenated_data_full',
                         inputs = 'data/03_concatenated_data/full_mouse/',
                         out    = 'data/04_integrated_data/full_mouse/',
                         scores = 'data/05_evaluated_score/full_mouse/'),
  dissected_mouse = list(concat = 'data/03_concatenated_data/03_concatenated_data_subset',
                         inputs = 'data/03_concatenated_data/dissected_mouse/',
                         out    = 'data/04_integrated_data/dissected_mouse/',
                         scores = 'data/05_evaluated_score/dissected_mouse/'),
  genes2          = list(concat = 'data/genes2/03_concatenated_data',
                         inputs = 'data/genes2/',
                         out    = 'data/genes2/',
                         scores = 'data/genes2/'),
  genes4          = list(concat = 'data/genes4/03_concatenated_data',
                         inputs = 'data/genes4/',
                         out    = 'data/genes4/',
                         scores = 'data/genes4/'),
  genes6          = list(concat = 'data/genes6/03_concatenated_data',
                         inputs = 'data/genes6/',
                         out    = 'data/genes6/',
                         scores = 'data/genes6/'),
  genes8          = list(concat = 'data/genes8/03_concatenated_data',
                         inputs = 'data/genes8/',
                         out    = 'data/genes8/',
                         scores = 'data/genes8/')
)
for (d in datasets) for (p in c(d$out, d$scores)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

# dataset names on the command line run only those, e.g. Rscript 04.01_integrate_fastMNN.R genes2 genes4
run_only = commandArgs(trailingOnly = TRUE)
if (length(run_only) > 0) {
  stopifnot(all(run_only %in% names(datasets)))
  datasets = datasets[run_only]
}
# threads come from N_THREADS through default_threads() in source.R
