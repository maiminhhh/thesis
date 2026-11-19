#!/usr/bin/env bash
# After the mouse EBCC (section 3b) has finished: dissection record with EBCC labels (section 5b), then
# 03 -> 04 -> 05 -> figures (run_pipeline.sh). Do not edit while it runs.
set -u
cd "$(dirname "$0")/.."
until grep -qE "3b FINISHED|STOP" logs/02/mouse_3b.log 2>/dev/null; do sleep 300; done
if grep -q STOP logs/02/mouse_3b.log; then echo "STOP: mouse 3b failed, nothing started ($(date))"; exit 1; fi
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 N_THREADS=4
echo "== dissection map + record with EBCC labels ($(date))"
nice -n 10 Rscript pipeline/mouse_dissection_map.R ebcc_label ebcc data/02_annotated_data/mouse_dissection_candidates_final.csv > logs/02/mouse_dissection_map_ebcc.log 2>&1 || echo "note: EBCC map failed (see log)"
nice -n 10 Rscript pipeline/mouse_dissection_final.R > logs/02/mouse_dissection_final.log 2>&1 || { echo "STOP: dissection record failed"; exit 1; }
tail -14 logs/02/mouse_dissection_final.log
echo "== starting run_pipeline.sh ($(date))"
SKIP_WAIT=1 bash run_pipeline.sh > logs/pipeline.log 2>&1
