#!/usr/bin/env bash
# Run every 04 integration with one process per dataset, MAX_JOBS at a time.
#   bash run_04.sh                    # all datasets
#   MAX_JOBS=4 bash run_04.sh         # fewer parallel processes (an R job on full_mouse can use 20-40 GB)
#   DATASETS="genes8" bash run_04.sh  # one dataset only
# Stage 1: methods that only need 03 (+ LIGER UINMF, which makes the H matrix)
# Stage 2: the UINMF-H methods, which need stage 1's uinmf_H_raw.* files
# GPU jobs run alongside on their own GPU queue each (one model per GPU at a time):
#   GPU 0: SAMap on full_mouse
#   GPU 1: scVI then scANVI on every dataset, then SAMap on the dissected-based datasets
# Each job gets THREADS (default 8) BLAS/OpenMP threads and runs under nice, so the server stays usable.
# Logs: logs/04/<script>_<dataset>.log
set -u
cd "$(dirname "$0")"

MAX_JOBS=${MAX_JOBS:-6}
# largest datasets first, so they don't become the long tail
DATASETS=${DATASETS:-"full_mouse dissected_mouse genes2 genes4 genes6 genes8"}   # as in integration_sets.R / datasets.py
THREADS=${THREADS:-8}   # per job; the server is shared, so 6 jobs x 8 threads by default
# N_THREADS must equal the OpenMP/BLAS cap (default_threads() in source.R reads it; Harmony refuses a mismatch)
export OPENBLAS_NUM_THREADS=$THREADS OMP_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS N_THREADS=$THREADS N_CPU=$THREADS
export RETICULATE_PYTHON=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python
PY_TORCH=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python
PY_SAMAP=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/samap/bin/python
LOG=logs/04
mkdir -p "$LOG"

# run "<interpreter> <script>" for every dataset, MAX_JOBS processes at a time
run_stage() {
  for g in $DATASETS; do
    for cmd in "$@"; do
      script=${cmd##* }
      echo "nice -n 10 $cmd $g > $LOG/${script%.*}_$g.log 2>&1"
    done
  done | xargs -P "$MAX_JOBS" -I{} bash -c '{} && echo "done:   {}" || echo "FAILED: {}"'
}

# split the datasets for the two GPU queues
FULL_DS=$(for g in $DATASETS; do case $g in full_*) echo -n "$g ";; esac; done)
DISS_DS=$(for g in $DATASETS; do case $g in full_*) ;; *) echo -n "$g ";; esac; done)

gpu_job() {  # gpu_job <gpu> <python> <script> <dataset>: one process per dataset
  local gpu=$1 py=$2 script=$3 g=$4
  CUDA_VISIBLE_DEVICES=$gpu nice -n 10 $py $script $g > "$LOG/${script%.*}_$g.log" 2>&1 \
    && echo "done:   $script $g" || { echo "FAILED: $script $g"; return 1; }
}

( for g in $FULL_DS; do gpu_job 0 $PY_SAMAP 04.02_integrate_SAMap.py $g; done ) &
( for g in $DATASETS; do   # scANVI reuses the scVI model, so it runs right after scVI per dataset
    gpu_job 1 $PY_TORCH 04.07_integrate_scVI.py $g && gpu_job 1 $PY_TORCH 04.08_integrate_scANVI.py $g
  done
  for g in $DISS_DS; do gpu_job 1 $PY_SAMAP 04.02_integrate_SAMap.py $g; done ) &

echo "== 04 stage 1 ($(date))"
run_stage "Rscript 04.11_integrate_LIGER_UINMF.R" \
          "Rscript 04.01_integrate_fastMNN.R" "Rscript 04.03_integrate_CCA.R" \
          "Rscript 04.04_integrate_Harmony.R" "Rscript 04.06_integrate_RPCA.R" \
          "Rscript 04.10_integrate_LIGER.R" \
          "$PY_TORCH 04.05_integrate_Scanorama.py" "$PY_TORCH 04.09_integrate_BBKNN.py"

echo "== 04 stage 2: UINMF H ($(date))"
run_stage "Rscript 04.12_integrate_uinmf_fastMNN.R" "Rscript 04.13_integrate_uinmf_CCA.R" \
          "Rscript 04.14_integrate_uinmf_harmony.R" "Rscript 04.16_integrate_uinmf_RPCA.R" \
          "$PY_TORCH 04.15_integrate_uinmf_scanorama.py" "$PY_TORCH 04.17_integrate_uinmf_bbknn.py"

wait   # GPU queues
echo "== all 04 jobs finished ($(date)); check $LOG for FAILED lines above"
