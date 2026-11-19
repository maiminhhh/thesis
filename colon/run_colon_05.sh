#!/usr/bin/env bash
# The tail of run_colon_frac.sh (all of 05, then figures), restarted on 2026-10-24 with more, thinner jobs:
# the scib scoring jobs run on ~1 core each, so 12 jobs x 2 threads keeps the 24-thread cap and runs ~2x faster.
# Runs in screen `colon_frac` (same name, so the label-check and dissection waiters keep waiting for it).
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
STATUS=logs/pipeline_status.txt
stage() { echo "$(date '+%F %T') $*" | tee -a $STATUS; }
trap 'stage "05 runner killed by a signal"; exit 1' TERM INT HUP

stage "05 restarted with 12 jobs x 2 threads (all methods, all datasets)"
mkdir -p logs/05; : > logs/05/status.txt
N_JOBS=12 N_THREADS=2 bash run_05.sh > logs/run_05.log 2>&1
stage "05 done: $(grep -vc 'exit=0' logs/05/status.txt) failed jobs (logs/05/status.txt)"

# figures (shared ../analysis/visualisations.R, colon only): fresh UMAP and spatial folders, since it skips existing files
rm -rf ../figures/colon/umap ../figures/colon/spatial
(cd .. && OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 N_THREADS=8 Rscript analysis/visualisations.R colon) > logs/visualisations.log 2>&1 ||
  { stage "FAILED figures (logs/visualisations.log)"; exit 1; }
stage "figures done (figures/colon/)"
