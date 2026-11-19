#!/usr/bin/env bash
# Score every 04 integration with scib, one job per (script, dataset), several at a time.
# A failed job doesn't stop the others; each job logs to logs/05/<script>__<dataset>.log
#
# usage: bash run_05.sh [dataset ...]        (default: every dataset in integration_datasets.py)
#   N_JOBS=8 N_THREADS=4 bash run_05.sh genes2 genes4
#   DRY_RUN=1 bash run_05.sh     (print the jobs instead of running them)
#   LABEL_SET=best_annotator LABEL_COLS='hsapiens=celltype_seurat,mmusculus=ebcc_label' bash run_05.sh
#        (labels from other 03 columns; scores in label_<set>/, logs in logs/05_<set>/)
#   LOGS_DIR=logs/05_x bash run_05.sh   (log folder override)
#   JOB_LIST=jobs.txt bash run_05.sh    (only these "script dataset" pairs)
#   SCRIPTS="05.01_evaluate_fastMNN.py 05.03_evaluate_CCA.py" bash run_05.sh   (only these scripts)
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/fetal_liver

N_JOBS=${N_JOBS:-8}          # jobs at once
N_THREADS=${N_THREADS:-4}    # threads per job
PYTHON=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/scib/bin/python
export OMP_NUM_THREADS=$N_THREADS OPENBLAS_NUM_THREADS=$N_THREADS MKL_NUM_THREADS=$N_THREADS \
       NUMEXPR_NUM_THREADS=$N_THREADS NUMBA_NUM_THREADS=$N_THREADS

DATASETS=("$@")
[ ${#DATASETS[@]} -eq 0 ] && DATASETS=(full_mouse dissected_mouse genes2 genes4 genes6 genes8)

SCRIPTS=${SCRIPTS:-$(ls 05.*.py)}

# wait for any 04 run of this project, so no method is scored before its 04 output is final
fetal_04_running() {
  for p in $(pgrep -f 'run_04.sh'); do
    [ "$(readlink /proc/$p/cwd 2>/dev/null)" == "$PWD" ] && return 0
  done
  return 1
}
[ -z "${DRY_RUN:-}" ] && while fetal_04_running; do echo "$(date +%T) waiting for run_04.sh"; sleep 120; done

LABEL_SET=${LABEL_SET:-ebcc}
LOGS=logs/05; [ "$LABEL_SET" != ebcc ] && LOGS=logs/05_$LABEL_SET
LOGS=${LOGS_DIR:-$LOGS}   # LOGS_DIR overrides the log folder
# one run per label set: if a run for this label set is already going (e.g. started early,
# alongside another run), wait for it and report its status instead of scoring everything twice
lock=logs/05_$LABEL_SET.lock
if [ -z "${DRY_RUN:-}" ] && [ -s $lock ] && kill -0 "$(head -1 $lock)" 2>/dev/null; then
  other_logs=$(sed -n 2p $lock)
  echo "a run for LABEL_SET=$LABEL_SET is already going (pid $(head -1 $lock), logs in $other_logs): waiting for it"
  while kill -0 "$(head -1 $lock)" 2>/dev/null; do sleep 120; done
  mkdir -p $LOGS
  [ "$other_logs" != "$LOGS" ] && cp "$other_logs/status.txt" $LOGS/status.txt
  echo "== summary"; sort $LOGS/status.txt
  echo "failed jobs: $(grep -vc 'exit=0' $LOGS/status.txt) (logs in $other_logs/)"
  exit 0
fi
[ -z "${DRY_RUN:-}" ] && printf '%s\n%s\n' $$ $LOGS > $lock
mkdir -p $LOGS
touch $LOGS/status.txt   # appended to, so separate batches share one summary

run_job() {
  local script=$1 dataset=$2 log="$LOGS/${1%.*}__$2.log" start=$SECONDS
  local cmd=("$PYTHON" -u "$script" "$dataset")   # -u: progress reaches the log as it happens
  [ -n "${DRY_RUN:-}" ] && cmd=(echo "${cmd[@]}")
  "${cmd[@]}" > "$log" 2>&1
  local rc=$?
  printf '%-40s %-16s exit=%d  %5d min\n' "$script" "$dataset" $rc $(( (SECONDS - start) / 60 )) | tee -a $LOGS/status.txt
}
export -f run_job
export LOGS PYTHON LABEL_SET DRY_RUN=${DRY_RUN:-}

# JOB_LIST: a file of "script dataset" lines, run instead of every SCRIPTS x dataset pair
{ if [ -n "${JOB_LIST:-}" ]; then cat "$JOB_LIST"
  else for s in $SCRIPTS; do for d in "${DATASETS[@]}"; do echo "$s $d"; done; done; fi; } |
  xargs -P "$N_JOBS" -L 1 bash -c 'run_job "$0" "$1"'

echo "== summary"
sort $LOGS/status.txt
echo "failed jobs: $(grep -vc 'exit=0' $LOGS/status.txt) (logs in $LOGS/)"
