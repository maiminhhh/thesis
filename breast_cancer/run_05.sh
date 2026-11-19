#!/bin/bash
# Run all 05 evaluations in parallel (needs every 04 output and eval_cells.txt).
#   ./run_05.sh                                   # all datasets, primary EBCC labels
#   JOBS=6 THREADS=8 ./run_05.sh
#   DATASETS="genes2 genes4" ./run_05.sh          # some datasets only
#   SCRIPTS="05.1[89]*.py 05.2*.py" ./run_05.sh  # only these scripts (e.g. the scVI-latent methods)
#   LABEL_SET=best_annotator LABEL_COLS='hsapiens=celltype_seurat,mmusculus=ebcc_label' JOBS=6 THREADS=8 ./run_05.sh
#                                                 # label sensitivity check (../EBCC_PIPELINE.md 9b)
# Scores: data/<dataset>/05.*.csv, or data/<dataset>/label_<set>/ for a LABEL_SET run.
# Logs: logs/05/ (logs/05_<set>/ for a LABEL_SET run)<script>_<dataset>.log; failures in <log dir>/FAILED.
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/breast_cancer
JOBS=${JOBS:-6}
THREADS=${THREADS:-8}
DATASETS=${DATASETS:-"full_mouse dissected_mouse genes2 genes4 genes6 genes8"}   # as in datasets.py
export PY_SCIB=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/scib/bin/python
export LABEL_SET=${LABEL_SET:-ebcc}
if [ "$LABEL_SET" = ebcc ]; then export LOG=logs/05; else export LOG=logs/05_$LABEL_SET LABEL_COLS; fi
mkdir -p "$LOG"
rm -f "$LOG/FAILED"
export OMP_NUM_THREADS=$THREADS OPENBLAS_NUM_THREADS=$THREADS MKL_NUM_THREADS=$THREADS NUMEXPR_NUM_THREADS=$THREADS N_CPU=$THREADS

for g in $DATASETS; do for f in $(ls ${SCRIPTS:-05.*.py}); do echo "$f $g"; done; done | xargs -P "$JOBS" -L1 bash -c '
  f=$0; g=$1; name=${f%.py}_$g
  echo "$(date +%T) start $name"
  nice -n 10 $PY_SCIB "$f" "$g" > "$LOG/$name.log" 2>&1
  rc=$?
  echo "$(date +%T) done  $name (exit $rc)"
  [ $rc -eq 0 ] || echo "$name" >> "$LOG/FAILED"'

if [ -s "$LOG/FAILED" ]; then echo "FAILED: $(tr '\n' ' ' < "$LOG/FAILED")"; exit 1; fi
echo "ALL DONE ($LABEL_SET)"
