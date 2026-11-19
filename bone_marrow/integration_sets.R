# datasets every 04 integration script runs on: concatenated input, the
# per-species inputs (LIGER UINMF) and the output folder. Mirrors datasets.py.
integration_sets = list(
  full = list(
    concat = 'data/03_concatenated_data/03_concatenated_data_full.rds',
    human  = 'data/02_annotated_data/02_human_annotated.rds',
    mouse  = 'data/02_annotated_data/02_mouse_annotated.rds',
    out    = 'data/04_integrated_data/full_mouse/'),
  dissected = list(
    concat = 'data/03_concatenated_data/03_concatenated_data_subset.rds',
    human  = 'data/02_annotated_data/02_human_annotated.rds',
    mouse  = 'data/02_annotated_data/02_mouse_subset_annotated.rds',
    out    = 'data/04_integrated_data/dissected_mouse/')
)
for (i in c(2, 4, 6, 8)) {
  integration_sets[[paste0('genes', i)]] = list(
    concat = str_glue('data/genes{i}/03_concatenated_data.rds'),
    human  = str_glue('data/genes{i}/02_human_annotated.rds'),
    mouse  = str_glue('data/genes{i}/02_mouse_annotated.rds'),
    out    = str_glue('data/genes{i}/'))
}

# datasets to run this time: names given on the command line, or all of them
# (e.g. `Rscript 04.01_integrate_fastMNN.R genes2 genes4`)
run_sets = commandArgs(trailingOnly = TRUE)
if (length(run_sets) == 0) run_sets = names(integration_sets)
stopifnot(all(run_sets %in% names(integration_sets)))
