#!/bin/bash
# Run 03 and 04 in parallel, respecting their dependencies.
#   ./run_03_04.sh          # 03 then 04
#   ./run_03_04.sh 03       # only 03
#   ./run_03_04.sh 04       # only 04 (03 outputs must exist)
#   ./run_03_04.sh 04scvi   # only the scVI-latent methods (04.18-04.24)
# Tune for how busy the (shared) machine is:
#   JOBS=6 THREADS=8 ./run_03_04.sh
# DRY=1 prints the jobs without running them.
# Logs: logs/03_04/<job>.log; failed jobs are listed in logs/03_04/FAILED.
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/bone_marrow

PART=${1:-all}
JOBS=${JOBS:-6}          # concurrent CPU jobs
THREADS=${THREADS:-8}    # BLAS/OpenMP threads per job
PY_TORCH=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python   # scVI, scANVI, Scanorama, BBKNN
PY_SAMAP=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/samap/bin/python
DATASETS=${DATASETS:-"full dissected genes2 genes4 genes6 genes8"}   # override: DATASETS="dissected genes2" ./run_03_04.sh 04

export LOG=logs/03_04
mkdir -p "$LOG"
rm -f "$LOG/FAILED"
# N_THREADS feeds default_threads() in source.R (Harmony, LIGER, BiocParallel, UMAP);
# keep it equal to the OpenMP cap or Harmony refuses to start
export N_THREADS=$THREADS OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS NUMEXPR_NUM_THREADS=$THREADS

# job NAME CMD...: run CMD, log to $LOG/NAME.log, note failures
job() {
  local name=$1; shift
  if [ -n "${DRY:-}" ]; then echo "[dry] $name: CUDA=${CUDA_VISIBLE_DEVICES-unset} $*"; return 0; fi
  echo "$(date +%T) start $name"
  "$@" > "$LOG/$name.log" 2>&1
  local rc=$?
  echo "$(date +%T) done  $name (exit $rc)"
  [ $rc -eq 0 ] || echo "$name" >> "$LOG/FAILED"
  return $rc
}
export -f job

# pool N: run "NAME CMD..." lines from stdin, N at a time
pool() { xargs -r -P "$1" -L 1 bash -c 'job "$@"' _; }

# ---- 03 ----
if [ "$PART" = all ] || [ "$PART" = 03 ]; then
  job 03_inputs Rscript 03_concatenate_data.R inputs &
  job 03_full   Rscript 03_concatenate_data.R full &
  job 03_subset Rscript 03_concatenate_data.R subset   # the gene levels need this
  for p in 0.2 0.4 0.6 0.8; do echo "03_genes$p Rscript 03_concatenate_data.R $p"; done | pool 4
  wait
  if [ -s "$LOG/FAILED" ]; then echo "03 failed: $(cat "$LOG/FAILED"); not starting 04"; exit 1; fi
fi

# ---- 04 ----
if [ "$PART" = all ] || [ "$PART" = 04 ]; then
  # GPU lane 0: scVI then scANVI (scANVI loads the scVI model)
  ( for s in $DATASETS; do
      CUDA_VISIBLE_DEVICES=0 job 04.07_$s $PY_TORCH 04.07_integrate_scVI.py $s &&
      CUDA_VISIBLE_DEVICES=0 job 04.08_$s $PY_TORCH 04.08_integrate_scANVI.py $s
      # scVI latent without the species key: the input of 04.19-04.24
      CUDA_VISIBLE_DEVICES=0 job 04.18_$s $PY_TORCH 04.18_scVI_latent.py $s
    done ) &
  GPU0=$!
  # GPU lane 1: SAMap
  ( for s in $DATASETS; do
      CUDA_VISIBLE_DEVICES=1 job 04.02_$s $PY_SAMAP 04.02_integrate_SAMap.py $s
    done ) &

  # CPU phase 1: LIGER UINMF (H matrix + eval_cells) first, then the gene-based methods
  {
    for s in $DATASETS; do echo "04.11_$s Rscript 04.11_integrate_LIGER_UINMF.R $s"; done
    for s in $DATASETS; do
      for f in 04.01_integrate_fastMNN 04.03_integrate_CCA 04.04_integrate_Harmony \
               04.06_integrate_RPCA 04.10_integrate_LIGER; do
        # spectral UMAP init hangs on RPCA genes2: random start there (see 04.06)
        pre=""; [ "$f" = 04.06_integrate_RPCA ] && [ "$s" = genes2 ] && pre="env UMAP_INIT=random"
        echo "${f%%_*}_$s $pre Rscript $f.R $s"
      done
      echo "04.05_$s $PY_TORCH 04.05_integrate_Scanorama.py $s"
      echo "04.09_$s $PY_TORCH 04.09_integrate_BBKNN.py $s"
    done
  } | CUDA_VISIBLE_DEVICES="" pool "$JOBS"

  # CPU phase 2: methods built on the UINMF H matrix
  {
    for s in $DATASETS; do
      for f in 04.12_integrate_fastMNN_uinmf 04.13_integrate_CCA_uinmf \
               04.14_integrate_Harmony_uinmf 04.16_integrate_RPCA_uinmf; do
        echo "${f%%_*}_$s Rscript $f.R $s"
      done
      echo "04.15_$s $PY_TORCH 04.15_integrate_Scanorama_uinmf.py $s"
      echo "04.17_$s $PY_TORCH 04.17_integrate_BBKNN_uinmf.py $s"
    done
  } | CUDA_VISIBLE_DEVICES="" pool "$JOBS"

  # CPU phase 3: the same methods on the scVI latent (needs 04.18 from GPU lane 0)
  wait $GPU0
  {
    for s in $DATASETS; do
      for f in 04.19_integrate_fastMNN_scvi 04.20_integrate_CCA_scvi \
               04.21_integrate_Harmony_scvi 04.23_integrate_RPCA_scvi; do
        echo "${f%%_*}_$s Rscript $f.R $s"
      done
      echo "04.22_$s $PY_TORCH 04.22_integrate_Scanorama_scvi.py $s"
      echo "04.24_$s $PY_TORCH 04.24_integrate_BBKNN_scvi.py $s"
    done
  } | CUDA_VISIBLE_DEVICES="" pool "$JOBS"

  wait
fi

# ---- only the scVI-latent methods (04.18 on the GPU, then 04.19-04.24) ----
if [ "$PART" = 04scvi ]; then
  for s in $DATASETS; do echo "04.18_$s $PY_TORCH 04.18_scVI_latent.py $s"; done |
    CUDA_VISIBLE_DEVICES=0 pool 1
  {
    for s in $DATASETS; do
      [ -s "data/04_integrated_data/${s}_mouse/scvi_latent.csv" ] || [ -s "data/$s/scvi_latent.csv" ] || continue
      for f in 04.19_integrate_fastMNN_scvi 04.20_integrate_CCA_scvi \
               04.21_integrate_Harmony_scvi 04.23_integrate_RPCA_scvi; do
        echo "${f%%_*}_$s Rscript $f.R $s"
      done
      echo "04.22_$s $PY_TORCH 04.22_integrate_Scanorama_scvi.py $s"
      echo "04.24_$s $PY_TORCH 04.24_integrate_BBKNN_scvi.py $s"
    done
  } | CUDA_VISIBLE_DEVICES="" pool "$JOBS"
fi

if [ -s "$LOG/FAILED" ]; then echo "FAILED: $(tr '\n' ' ' < "$LOG/FAILED")"; exit 1; fi
echo "ALL DONE"
