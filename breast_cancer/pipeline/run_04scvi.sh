#!/usr/bin/env bash
# scVI-latent methods (../EBCC_PIPELINE.md 6b) for the datasets given (default: all but genes8, the smoke test):
# 04.18 per dataset on two GPU queues, and as each latent appears its 04.19-04.24 jobs on CPU (at most JOBS at once,
# THREADS each). Results: logs/04scvi/run.log ("<step> <dataset>: exit N", then "== all 04scvi finished").
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/breast_cancer
JOBS=${JOBS:-6}; T=${THREADS:-4}
GPU0=${GPU0:-"dissected_mouse genes2 full_mouse"}; GPU1=${GPU1:-"genes4 genes6"}
export RETICULATE_PYTHON=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python
export N_THREADS=$T N_CPU=$T OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T NUMEXPR_NUM_THREADS=$T
PY_TORCH=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python
L=logs/04scvi; R=$L/run.log; mkdir -p $L
job() { local step=$1 g=$2; shift 2; nice -n 10 "$@" > $L/${step}_$g.log 2>&1; echo "$(date '+%F %T') $step $g: exit $?" >> $R; }
gpu_queue() { local gpu=$1; shift; for g in "$@"; do CUDA_VISIBLE_DEVICES=$gpu job 04.18 $g $PY_TORCH 04.18_scVI_latent.py $g; done; }
echo "$(date '+%F %T') == start (GPU0: $GPU0 | GPU1: $GPU1)" >> $R
gpu_queue 0 $GPU0 & gpu_queue 1 $GPU1 &
todo="$GPU0 $GPU1"
while [ -n "$todo" ]; do
  left=""
  for g in $todo; do
    if grep -q "04.18 $g: exit 0" $R; then
      for c in "04.19 Rscript 04.19_integrate_fastMNN_scvi.R" "04.20 Rscript 04.20_integrate_CCA_scvi.R" \
               "04.21 Rscript 04.21_integrate_Harmony_scvi.R" "04.23 Rscript 04.23_integrate_RPCA_scvi.R" \
               "04.22 $PY_TORCH 04.22_integrate_Scanorama_scvi.py" "04.24 $PY_TORCH 04.24_integrate_BBKNN_scvi.py"; do
        while [ $(jobs -rp | wc -l) -ge $((JOBS + 2)) ]; do sleep 30; done   # +2: the GPU queues
        set -- $c; CUDA_VISIBLE_DEVICES="" job $1 $g $2 $3 $g &
      done
    elif grep -qE "04.18 $g: exit [1-9]" $R; then echo "$(date '+%F %T') SKIP: $g has no latent" >> $R
    else left="$left $g"; fi
  done
  todo=$(echo $left); [ -n "$todo" ] && sleep 60
done
wait
echo "$(date '+%F %T') == all 04scvi finished" >> $R
