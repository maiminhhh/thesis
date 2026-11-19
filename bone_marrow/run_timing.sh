#!/bin/bash
# Controlled runtime benchmark: one method at a time, integration step only,
# N_THREADS CPU threads, GPU methods on one GPU; the order is shuffled per repeat
# (seeded) so changes in server load spread over methods.
#   ./run_timing.sh                       # dissected + full, 3 repeats
#   DATASETS=genes8 REPS=1 ./run_timing.sh
# Results: data/08_timing/timing.csv; logs: logs/08_timing/.
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/bone_marrow
DATASETS=${DATASETS:-"dissected full"}
REPS=${REPS:-3}
export N_THREADS=${N_THREADS:-8}
export OMP_NUM_THREADS=$N_THREADS OPENBLAS_NUM_THREADS=$N_THREADS MKL_NUM_THREADS=$N_THREADS NUMEXPR_NUM_THREADS=$N_THREADS
export RETICULATE_PYTHON=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python
PY_TORCH=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python
PY_SAMAP=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/samap/bin/python
LOG=logs/08_timing; mkdir -p "$LOG"
L=$LOG/run.log; : > $L
step() { echo "===== $1 $(date)" >> $L; }

R_METHODS=("FastMNN" "CCA" "Harmony" "RPCA" "LIGER" "LIGER UINMF"
           "FastMNN UINMF" "CCA UINMF" "Harmony UINMF" "RPCA UINMF"
           "FastMNN scVI-latent" "CCA scVI-latent" "Harmony scVI-latent" "RPCA scVI-latent")
PY_METHODS=("Scanorama" "BBKNN" "Scanorama UINMF" "BBKNN UINMF"
            "Scanorama scVI-latent" "BBKNN scVI-latent" "scVI" "scANVI" "scVI latent" "SAMap")

n_fail=0
for rep in $(seq 1 "$REPS"); do
  jobs=()
  for s in $DATASETS; do
    for m in "${R_METHODS[@]}" "${PY_METHODS[@]}"; do jobs+=("$s|$m"); done
  done
  step "repeat $rep: ${#jobs[@]} runs"
  printf '%s\n' "${jobs[@]}" | shuf --random-source=<(yes $rep) | while IFS='|' read -r s m; do
    tag="${m// /_}_${s}_r$rep"
    if printf '%s\n' "${R_METHODS[@]}" | grep -qx "$m"; then
      cmd=(Rscript 08_timing.R "$m" "$s" "$rep")
    elif [ "$m" = SAMap ]; then
      cmd=(env CUDA_VISIBLE_DEVICES=0 $PY_SAMAP 08_timing.py "$m" "$s" "$rep")
    else
      cmd=(env CUDA_VISIBLE_DEVICES=0 $PY_TORCH 08_timing.py "$m" "$s" "$rep")
    fi
    echo "$(date +%T) start $tag" >> $L
    "${cmd[@]}" > "$LOG/$tag.log" 2>&1
    rc=$?
    echo "$(date +%T) done  $tag (exit $rc)" >> $L
    [ $rc -eq 0 ] || echo "$tag" >> "$LOG/FAILED"
  done
done
if [ -s "$LOG/FAILED" ]; then step "FAILED: $(wc -l < "$LOG/FAILED") runs (see $LOG/FAILED)"; exit 1; fi
step "ALL DONE"
