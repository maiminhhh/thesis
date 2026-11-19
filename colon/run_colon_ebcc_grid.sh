#!/usr/bin/env bash
# EBCC procedure (../EBCC_PIPELINE.md) for colon, then the 03/04/05 pipeline and figures:
#   1. rerun the 4 annotators that use the human (CELLxGENE) reference (expm1(X) fix)
#   2. (M, frac) held-out grid, one job per (species, frac), each over all M
#   3. per species: pick (M, frac), write it into the 02 script, final fit (labels changed, so always)
#   4. run_colon_pipeline.sh: 02 h5ad export -> 03 -> 04 -> 05 -> figures
# Run it in screen:  screen -dmS colon_ebcc bash run_colon_ebcc_grid.sh
# Progress: logs/pipeline_status.txt; logs in logs/ebcc_grid/
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
G=logs/ebcc_grid
STATUS=logs/pipeline_status.txt
stage() { echo "$(date '+%F %T') $*" | tee -a $STATUS; }
trap 'stage "EBCC runner killed by a signal"; exit 1' TERM INT HUP
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1 N_THREADS=$1; }

stage "re-annotation started (human same-species, mouse cross-species)"
threads 8
for s in human mouse; do Rscript $G/reannotate.R $s > $G/reannotate_$s.log 2>&1 & done
fail=0; for j in $(jobs -p); do wait $j || fail=1; done
[ $fail = 1 ] && { stage "FAILED re-annotation ($G/reannotate_*.log)"; exit 1; }
stage "re-annotation done"

stage "EBCC grid started (10 jobs x 3 restarts)"
threads 1   # run_ebcc runs the 3 restarts as parallel processes
for s in human mouse; do for f in 0 0.25 0.5 1 2; do echo "$s $f"; done; done |
  xargs -P 10 -L 1 bash -c 'Rscript '$G'/ebcc_grid_part.R "$0" "$1" > '$G'/part_$0_$1.log 2>&1 || echo "FAILED $0 frac=$1"' \
  > $G/failed.txt
[ -s $G/failed.txt ] && { stage "FAILED EBCC grid: $(tr '\n' ' ' < $G/failed.txt)(logs in $G/)"; exit 1; }
stage "EBCC grid done"

for s in human mouse; do
  Rscript $G/ebcc_grid_choose.R $s > $G/choose_$s.log 2>&1 || { stage "FAILED choosing $s (M, frac) ($G/choose_$s.log)"; exit 1; }
  stage "$(grep 'chosen M' $G/choose_$s.log)"
done

touch $G/.fit_started
stage "final EBCC fits started (10 restarts each, in parallel)"
for s in human mouse; do Rscript $G/ebcc_grid_finish.R $s > $G/finish_$s.log 2>&1 & done
fail=0; for j in $(jobs -p); do wait $j || fail=1; done
for s in human mouse; do
  [ data/02_annotated_data/02_${s}_annotated.rds -nt $G/.fit_started ] || fail=1
done
[ $fail = 1 ] && { stage "FAILED final EBCC fit ($G/finish_*.log)"; exit 1; }
stage "final EBCC fits done"

exec bash run_colon_pipeline.sh
