#!/usr/bin/env bash
# After the 02 EBCC fits: export the 02 h5ad files, then run 03, all of 04 and all of 05.
# Each stage writes one line to logs/pipeline_status.txt; a failed stage stops the chain.
#
# usage: bash run_pipeline.sh [pid ...]   (wait for these processes, e.g. the EBCC fits, first)
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
mkdir -p logs
STATUS=logs/pipeline_status.txt
stage() { echo "$(date '+%F %T') $*" | tee -a $STATUS; }
trap 'stage "pipeline runner killed by a signal"; exit 1' TERM INT HUP

for pid in "$@"; do
  while kill -0 "$pid" 2>/dev/null; do sleep 120; done
done
stage "EBCC fits finished"

Rscript export_02_h5ad.R > logs/export_02_h5ad.log 2>&1 || { stage "FAILED 02 export/check (logs/export_02_h5ad.log)"; exit 1; }
stage "02 h5ad exported"

# old outputs out of the way, so a failed method can't leave a stale result for 05 to score
OLD=data/old_runs_$(date +%F_%H%M)
mkdir -p $OLD/04_integrated_data $OLD/05_evaluated_score
for p in data/03_concatenated_data.rds data/03_concatenated_data.h5ad \
         data/03_concatenated_data_subset.rds data/03_concatenated_data_subset.h5ad \
         data/03_concatenated_data data/genes2 data/genes4 data/genes6 data/genes8; do
  [ -e "$p" ] && mv "$p" $OLD/
done
for p in data/04_integrated_data/*; do [ -e "$p" ] && mv "$p" $OLD/04_integrated_data/; done
for p in data/05_evaluated_score/*; do [ -e "$p" ] && mv "$p" $OLD/05_evaluated_score/; done
stage "old outputs moved to $OLD"

Rscript -e 'knitr::purl("03_concatenate.qmd", output = "logs/03_concatenate.R", quiet = TRUE, documentation = 0)' &&
  N_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 Rscript logs/03_concatenate.R > logs/03_concatenate.log 2>&1 ||
  { stage "FAILED 03 (logs/03_concatenate.log)"; exit 1; }
stage "03 done"

N_JOBS=4 N_THREADS=8 bash run_04.sh > logs/run_04.log 2>&1
stage "04 done: $(grep -vc 'exit=0' logs/04/status.txt) failed jobs (logs/04/status.txt)"

N_JOBS=6 N_THREADS=4 bash run_05.sh > logs/run_05.log 2>&1
stage "05 done: $(grep -vc 'exit=0' logs/05/status.txt) failed jobs (logs/05/status.txt)"

# figures (shared ../analysis/visualisations.R, colon only): fresh UMAP and spatial folders, since it skips existing files
rm -rf ../figures/colon/umap ../figures/colon/spatial
(cd .. && OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 N_THREADS=8 Rscript analysis/visualisations.R colon) > logs/visualisations.log 2>&1 ||
  { stage "FAILED figures (logs/visualisations.log)"; exit 1; }
stage "figures done (figures/colon/)"
