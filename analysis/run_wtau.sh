#!/bin/bash
# Re-make the tissue-agreement panels with weighted Kendall tau (no rescoring): analysis/benchmark_figures.R (+MAP +SAS, and
# graph-metrics variant) and analysis/mcb_analysis.R. analysis/benchmark_figures.R is run from a copy with source('library.R') replaced by
# library(tidyverse) (Seurat no longer loads under R 4.5; nothing else from library.R is used).
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}
LOG=figures/wtau_run.log; R=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/r45/bin/Rscript
export R_LIBS_SITE= R_LIBS_USER=
: > $LOG
# The two copies were made by hand for the original run and not kept; they are rebuilt here from the description above
# (the graph-metrics variant = the graph_metrics_only switch of analysis/benchmark_figures.R set to TRUE).
TMP=$(mktemp -d)
sed "s/^source('library.R')/library(tidyverse)/" analysis/benchmark_figures.R > $TMP/bf_run.R
sed "s/^graph_metrics_only = FALSE/graph_metrics_only = TRUE /" $TMP/bf_run.R > $TMP/bf_run_graph.R
echo "===== analysis/benchmark_figures.R (+MAP +SAS)" >> $LOG; $R --vanilla $TMP/bf_run.R >> $LOG 2>&1 || { echo "===== FAILED: benchmark_figures" >> $LOG; exit 1; }
echo "===== analysis/benchmark_figures.R (graph metrics only)" >> $LOG; $R --vanilla $TMP/bf_run_graph.R >> $LOG 2>&1 || { echo "===== FAILED: benchmark_figures graph" >> $LOG; exit 1; }
echo "===== analysis/mcb_analysis.R" >> $LOG; $R --vanilla analysis/mcb_analysis.R >> $LOG 2>&1 || { echo "===== FAILED: mcb_analysis" >> $LOG; exit 1; }
echo "===== ALL DONE" >> $LOG
