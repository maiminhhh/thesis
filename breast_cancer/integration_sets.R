# datasets every 04 integration script runs on; each is a folder data/<name>/ holding
# 02_human_annotated, 02_mouse_annotated and 03_concatenated_data (.rds + .h5ad). Mirrors datasets.py.
integration_sets = c('full_mouse', 'dissected_mouse', 'genes2', 'genes4', 'genes6', 'genes8')

# datasets to run this time: names given on the command line, or all of them
# (e.g. `Rscript 04.01_integrate_fastMNN.R genes2 genes4`)
run_sets = commandArgs(trailingOnly = TRUE)
if (length(run_sets) == 0) run_sets = integration_sets
stopifnot(all(run_sets %in% integration_sets))
