#!/bin/bash
# Run all 05 evaluations in parallel (needs every 04 output and eval_cells.txt).
#   ./run_05.sh                 # all datasets
#   DATASETS="dissected genes2" ./run_05.sh   # only these
#   SCRIPTS="05.1[89]*.py 05.2*.py" ./run_05.sh  # only these scripts (e.g. the scVI-latent methods)
#   JOBS=6 THREADS=8 ./run_05.sh
# Logs: logs/05/<script>.log; failed scripts are listed in logs/05/FAILED.
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/bone_marrow
JOBS=${JOBS:-6}
THREADS=${THREADS:-8}
export DATASETS=${DATASETS:-}
export PY_SCIB=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/scib/bin/python
# LABEL_SET=<set> LABEL_COLS=... (see utils.py) scores against another label set; logs go to logs/05_<set>
export LOG=logs/05${LABEL_SET:+_$LABEL_SET}
mkdir -p "$LOG" data/05_evaluated_score/full_mouse data/05_evaluated_score/dissected_mouse
rm -f "$LOG/FAILED"
export OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS NUMEXPR_NUM_THREADS=$THREADS

ls ${SCRIPTS:-05.*.py} | xargs -P "$JOBS" -I{} bash -c '
  f={}; name=${f%.py}
  echo "$(date +%T) start $name"
  $PY_SCIB "$f" $DATASETS > "$LOG/$name.log" 2>&1
  rc=$?
  echo "$(date +%T) done  $name (exit $rc)"
  [ $rc -eq 0 ] || echo "$name" >> "$LOG/FAILED"'

if [ -s "$LOG/FAILED" ]; then echo "FAILED: $(tr '\n' ' ' < "$LOG/FAILED")"; exit 1; fi
echo "ALL DONE"
