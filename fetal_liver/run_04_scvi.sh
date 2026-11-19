#!/usr/bin/env bash
# The scVI-latent methods (../EBCC_PIPELINE.md section 6b): 04.18 scVI without the species key
# on the GPUs (datasets split across GPU 0 and 1), then 04.19-04.24 on that latent, several at
# a time. Logs and the status file are shared with run_04.sh (logs/04/).
# usage: bash run_04_scvi.sh [dataset ...]     (default: every dataset)
#   N_JOBS=4 N_THREADS=4 bash run_04_scvi.sh genes8
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/fetal_liver

N_JOBS=${N_JOBS:-4}
N_THREADS=${N_THREADS:-4}
CONDA=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}
export RETICULATE_PYTHON=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python
export N_THREADS OMP_NUM_THREADS=$N_THREADS OPENBLAS_NUM_THREADS=$N_THREADS \
       MKL_NUM_THREADS=$N_THREADS NUMEXPR_NUM_THREADS=$N_THREADS NUMBA_NUM_THREADS=$N_THREADS

DATASETS=("$@")
[ ${#DATASETS[@]} -eq 0 ] && DATASETS=(full_mouse dissected_mouse genes2 genes4 genes6 genes8)
LOGS=logs/04
mkdir -p $LOGS
echo "# run $(date '+%F %T'): scVI-latent ${DATASETS[*]}" >> $LOGS/status.txt

run_job() {
  local script=$1 dataset=$2 log="$LOGS/${1%.*}__$2.log" start=$SECONDS
  case $script in
    *.R)  cmd=(Rscript "$script" "$dataset") ;;
    *.py) cmd=("$CONDA/torch/bin/python" -u "$script" "$dataset") ;;
  esac
  timeout 8h "${cmd[@]}" > "$log" 2>&1
  local rc=$?
  printf '%-40s %-16s exit=%d  %5d min\n' "$script" "$dataset" $rc $(( (SECONDS - start) / 60 )) | tee -a $LOGS/status.txt
}
export -f run_job
export LOGS CONDA

# 04.18 on the GPUs: alternate datasets between GPU 0 and 1, one at a time per GPU
lane() { local gpu=$1; shift; for d in "$@"; do CUDA_VISIBLE_DEVICES=$gpu run_job 04.18_scVI_latent.py $d; done; }
lane0=(); lane1=()
for i in "${!DATASETS[@]}"; do
  if (( i % 2 == 0 )); then lane0+=("${DATASETS[$i]}"); else lane1+=("${DATASETS[$i]}"); fi
done
echo "== 04.18 scVI latent on the GPUs"
lane 0 "${lane0[@]}" & gpu0=$!
[ ${#lane1[@]} -gt 0 ] && { lane 1 "${lane1[@]}" & gpu1=$!; }
wait $gpu0; [ -n "${gpu1:-}" ] && wait $gpu1

echo "== 04.19-04.24 on the scVI latent, $N_JOBS at a time"
for s in 04.19_integrate_fastMNN_scvi.R 04.20_integrate_CCA_scvi.R 04.21_integrate_Harmony_scvi.R \
         04.22_integrate_Scanorama_scvi.py 04.23_integrate_RPCA_scvi.R 04.24_integrate_BBKNN_scvi.py; do
  for d in "${DATASETS[@]}"; do
    out=$(case $d in full_mouse|dissected_mouse) echo data/04_integrated_data/$d/ ;; *) echo data/$d/ ;; esac)
    [ -s "${out}scvi_latent.csv" ] && echo "$s $d" || echo "skipping $s $d: no ${out}scvi_latent.csv" >&2
  done
done | CUDA_VISIBLE_DEVICES="" xargs -r -P "$N_JOBS" -L 1 bash -c 'run_job "$0" "$1"'

echo "== summary (this run)"
this_run=$(tail -n "+$(grep -n '^# run ' $LOGS/status.txt | tail -1 | cut -d: -f1)" $LOGS/status.txt | grep -v '^#')
echo "$this_run" | sort
echo "failed jobs: $(echo "$this_run" | grep -c 'exit=[1-9]')"
