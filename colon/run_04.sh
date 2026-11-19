#!/usr/bin/env bash
# Run every 04 integration, one job per (script, dataset), several at a time.
#   phase 1: methods that only need 03 outputs, plus 04.11 (LIGER UINMF)
#            GPU methods (scVI, scANVI, SAMap) run alongside on one lane per GPU
#   phase 2: methods that need 04.11's H matrix (04.12-04.17)
# A failed job doesn't stop the others; each job logs to logs/04/<script>__<dataset>.log
#
# usage: bash run_04.sh [dataset ...]        (default: every dataset in integration_datasets.R)
#   N_JOBS=6 N_THREADS=8 bash run_04.sh genes2 genes4
#   DRY_RUN=1 bash run_04.sh     (print the jobs instead of running them)
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon

N_JOBS=${N_JOBS:-6}              # CPU jobs at once
N_THREADS=${N_THREADS:-8}        # threads per job
JOB_TIMEOUT=${JOB_TIMEOUT:-12h}  # a job still running after this is killed (exit 124)
CONDA=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}
export N_THREADS OMP_NUM_THREADS=$N_THREADS OPENBLAS_NUM_THREADS=$N_THREADS \
       MKL_NUM_THREADS=$N_THREADS NUMEXPR_NUM_THREADS=$N_THREADS NUMBA_NUM_THREADS=$N_THREADS

DATASETS=("$@")
[ ${#DATASETS[@]} -eq 0 ] && DATASETS=(full_mouse dissected_mouse genes2 genes4 genes6 genes8)

LOGS=logs/04
mkdir -p $LOGS
: > $LOGS/status.txt

# run one script on one dataset with the right interpreter, recording its exit status
run_job() {
  local script=$1 dataset=$2 log="$LOGS/${1%.*}__$2.log" start=$SECONDS
  case $script in
    *.R)            cmd=(Rscript "$script" "$dataset") ;;
    04.02_*)        cmd=("$CONDA/samap/bin/python" "$script" "$dataset") ;;
    *.py)           cmd=("$CONDA/torch/bin/python" "$script" "$dataset") ;;
  esac
  [ -n "${DRY_RUN:-}" ] && cmd=(echo "${CUDA_VISIBLE_DEVICES:+GPU$CUDA_VISIBLE_DEVICES }${cmd[@]}")
  timeout "$JOB_TIMEOUT" "${cmd[@]}" > "$log" 2>&1
  local rc=$?
  printf '%-40s %-16s exit=%d  %5d min\n' "$script" "$dataset" $rc $(( (SECONDS - start) / 60 )) | tee -a $LOGS/status.txt
}
export -f run_job
export LOGS CONDA JOB_TIMEOUT DRY_RUN=${DRY_RUN:-}

# one GPU lane: its share of the datasets, GPU methods one after another
gpu_lane() {
  local gpu=$1; shift
  for d in "$@"; do
    for s in 04.07_integrate_scVI.py 04.08_integrate_scANVI.py 04.02_integrate_SAMap.py; do
      CUDA_VISIBLE_DEVICES=$gpu run_job $s $d
    done
  done
}

lane0=(); lane1=()
for i in "${!DATASETS[@]}"; do
  if (( i % 2 == 0 )); then lane0+=("${DATASETS[$i]}"); else lane1+=("${DATASETS[$i]}"); fi
done
gpu_lane 0 "${lane0[@]}" & gpu0=$!
[ ${#lane1[@]} -gt 0 ] && { gpu_lane 1 "${lane1[@]}" & gpu1=$!; }

# 04.11 first: the UINMF methods wait for it
phase1=(04.11_integrate_LIGER_UINMF.R 04.03_integrate_CCA.R 04.06_integrate_RPCA.R
        04.10_integrate_LIGER.R 04.01_integrate_fastMNN.R 04.04_integrate_Harmony.R
        04.05_integrate_Scanorama.py 04.09_integrate_BBKNN.py)
phase2=(04.13_integrate_uinmf_CCA.R 04.16_integrate_uinmf_RPCA.R 04.12_integrate_uinmf_fastMNN.R
        04.14_integrate_uinmf_harmony.R 04.15_integrate_uinmf_scanorama.py 04.17_integrate_uinmf_bbknn.py)

jobs_for() { for s in "$@"; do for d in "${DATASETS[@]}"; do echo "$s $d"; done; done; }

echo "== phase 1: ${#phase1[@]} methods x ${#DATASETS[@]} datasets, $N_JOBS at a time, $N_THREADS threads each"
jobs_for "${phase1[@]}" | xargs -P "$N_JOBS" -L 1 bash -c 'run_job "$0" "$1"'

echo "== phase 2: UINMF methods"
# a dataset whose 04.11 failed has no H matrix: skip its UINMF methods
ok=()
for d in "${DATASETS[@]}"; do
  if grep -q "^04.11_integrate_LIGER_UINMF.R *$d *exit=0" $LOGS/status.txt; then ok+=("$d")
  else echo "skipping UINMF methods for $d: 04.11 failed"; fi
done
[ ${#ok[@]} -gt 0 ] && DATASETS=("${ok[@]}") && \
  jobs_for "${phase2[@]}" | xargs -P "$N_JOBS" -L 1 bash -c 'run_job "$0" "$1"'

echo "== waiting for GPU lanes"
wait $gpu0; [ -n "${gpu1:-}" ] && wait $gpu1

echo "== summary"
sort $LOGS/status.txt
echo "failed jobs: $(grep -vc 'exit=0' $LOGS/status.txt) (logs in $LOGS/)"
