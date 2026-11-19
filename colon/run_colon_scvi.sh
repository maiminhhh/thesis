#!/usr/bin/env bash
# Methods on an scVI latent (../EBCC_PIPELINE.md 6b): 04.18 (scVI without the species key) now, on both GPUs,
# alongside the main 05; then, after run_colon_05.sh has drawn the figures, 04.19-04.24, their scores
# (05.18-05.24) and the figures again.
# Run in screen:  screen -dmS colon_scvi bash run_colon_scvi.sh
# Job logs: logs/04scvi/; job status lines: logs/04/status_scvi.txt (05: logs/05/status.txt), so the watchdog sees them.
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
STATUS=logs/pipeline_status.txt
LOG=logs/04scvi
PY=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python
mkdir -p $LOG
stage() { echo "$(date '+%F %T') $*" | tee -a $STATUS; }
trap 'stage "scVI-latent runner killed by a signal"; exit 1' TERM INT HUP
export RETICULATE_PYTHON=$PY LOG
export N_THREADS=4 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
start_line=$(( $(wc -l < $STATUS) + 1 ))

# job NAME CMD...: run CMD, log to $LOG/NAME.log, one status line in logs/04/status_scvi.txt
job() {
  local name=$1 start=$SECONDS; shift
  "$@" > "$LOG/$name.log" 2>&1
  local rc=$?
  printf '%-40s %-16s exit=%d  %5d min\n' "${name%%__*}" "${name##*__}" $rc $(( (SECONDS - start) / 60 )) >> logs/04/status_scvi.txt
  return $rc
}
export -f job

# ---- 04.18: scVI latent, datasets split over both GPUs (skip any already done) ----
stage "scVI-latent 04.18 started (both GPUs)"
lane() { local gpu=$1; shift
  for d in "$@"; do
    [ -s "$(Rscript -e "source('integration_datasets.R'); cat(datasets[['$d']]\$out)")04.18_scVI_latent.h5ad" ] && continue
    CUDA_VISIBLE_DEVICES=$gpu job "04.18_scVI_latent.py__$d" $PY 04.18_scVI_latent.py $d
  done; }
lane 0 full_mouse genes2 genes6 & lane 1 dissected_mouse genes4 genes8 & wait
stage "scVI-latent 04.18 done: $(grep '04.18' logs/04/status_scvi.txt | grep -vc 'exit=0') failed (logs/04/status_scvi.txt)"

# ---- wait for the main 05 run and its figures ----
until tail -n +$start_line $STATUS | grep -q 'figures done'; do
  screen -ls | grep -q colon_frac || { stage "scVI-latent runner: colon_frac exited without figures, stopping"; exit 1; }
  sleep 300
done

# ---- 04.19-04.24 on the latent (CPU), 4 jobs x 6 threads ----
export N_THREADS=6 OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 NUMEXPR_NUM_THREADS=6
stage "scVI-latent 04.19-04.24 started"
for d in full_mouse dissected_mouse genes2 genes4 genes6 genes8; do
  for f in 04.19_integrate_fastMNN_scvi.R 04.20_integrate_CCA_scvi.R 04.21_integrate_Harmony_scvi.R 04.23_integrate_RPCA_scvi.R; do
    echo "${f}__$d Rscript $f $d"; done
  echo "04.22_integrate_Scanorama_scvi.py__$d $PY 04.22_integrate_Scanorama_scvi.py $d"
  echo "04.24_integrate_BBKNN_scvi.py__$d $PY 04.24_integrate_BBKNN_scvi.py $d"
done | CUDA_VISIBLE_DEVICES="" xargs -r -P 4 -L 1 bash -c 'job "$@"' _
stage "scVI-latent 04.19-04.24 done: $(grep -E '04\.(19|2[0-4])' logs/04/status_scvi.txt | grep -vc 'exit=0') failed (logs/04/status_scvi.txt)"

# ---- 05.18-05.24 (12 jobs x 2 threads), then figures ----
stage "scVI-latent 05.18-05.24 started"
SCRIPTS="05.1[89]*.py 05.2*.py" N_JOBS=12 N_THREADS=2 bash run_05.sh > logs/run_05_scvi.log 2>&1
stage "scVI-latent 05 done: $(grep -E '05\.(1[89]|2[0-4])' logs/05/status.txt | grep -vc 'exit=0') failed (logs/05/status.txt)"
(cd .. && OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 N_THREADS=8 Rscript analysis/visualisations.R colon) > logs/visualisations_scvi.log 2>&1 ||
  { stage "FAILED figures after scVI-latent (logs/visualisations_scvi.log)"; exit 1; }
stage "figures redrawn with the scVI-latent methods (figures/colon/)"
