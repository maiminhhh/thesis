#!/usr/bin/env bash
# Score every 04 integration with scib, one job per (script, dataset), several at a time.
# A failed job doesn't stop the others; each job logs to logs/05/<script>__<dataset>.log
#
# usage: bash run_05.sh [dataset ...]        (default: every dataset in integration_datasets.py)
#   N_JOBS=8 N_THREADS=4 bash run_05.sh genes2 genes4
#   DRY_RUN=1 bash run_05.sh     (print the jobs instead of running them)
#   SCRIPTS="05.01_evaluate_fastMNN.py 05.03_evaluate_CCA.py" bash run_05.sh   (only these scripts)
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon

N_JOBS=${N_JOBS:-8}          # jobs at once
N_THREADS=${N_THREADS:-4}    # threads per job
PYTHON=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/scib/bin/python
export N_THREADS OMP_NUM_THREADS=$N_THREADS OPENBLAS_NUM_THREADS=$N_THREADS MKL_NUM_THREADS=$N_THREADS \
       NUMEXPR_NUM_THREADS=$N_THREADS NUMBA_NUM_THREADS=$N_THREADS

DATASETS=("$@")
[ ${#DATASETS[@]} -eq 0 ] && DATASETS=(full_mouse dissected_mouse genes2 genes4 genes6 genes8)

# 05.02_evaluate_scIB.py is an unfinished draft (two dictionaries, no scoring)
SCRIPTS=${SCRIPTS:-$(ls 05.0*.py 05.1*.py | grep -v scIB)}

# LABEL_SET / LABEL_COLS (utils.py): a non-default label set logs to logs/05_<set>/ (scores: label_<set>/ subfolders)
LABEL_SET=${LABEL_SET:-ebcc}
export LABEL_SET
LOGS=logs/05; [ "$LABEL_SET" != ebcc ] && LOGS=logs/05_$LABEL_SET
mkdir -p $LOGS
touch $LOGS/status.txt   # appended to, so separate batches share one summary

run_job() {
  local script=$1 dataset=$2 log="$LOGS/${1%.*}__$2.log" start=$SECONDS
  local cmd=("$PYTHON" "$script" "$dataset")
  [ -n "${DRY_RUN:-}" ] && cmd=(echo "${cmd[@]}")
  "${cmd[@]}" > "$log" 2>&1
  local rc=$?
  printf '%-40s %-16s exit=%d  %5d min\n' "$script" "$dataset" $rc $(( (SECONDS - start) / 60 )) | tee -a $LOGS/status.txt
}
export -f run_job
export LOGS PYTHON DRY_RUN=${DRY_RUN:-}

for s in $SCRIPTS; do for d in "${DATASETS[@]}"; do echo "$s $d"; done; done |
  xargs -P "$N_JOBS" -L 1 bash -c 'run_job "$0" "$1"'

echo "== summary"
sort $LOGS/status.txt
echo "failed jobs: $(grep -vc 'exit=0' $LOGS/status.txt) (logs in $LOGS/)"
