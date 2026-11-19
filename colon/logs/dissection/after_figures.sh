#!/usr/bin/env bash
# run the dissection report once run_colon_frac.sh has drawn the figures (it clears figures/colon/spatial first)
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
STATUS=logs/pipeline_status.txt
until grep -q 'figures done' <(tail -n +$1 $STATUS); do
  screen -ls | grep -q colon_frac || { echo "$(date '+%F %T') dissection report: colon_frac exited without figures, not run" >> $STATUS; exit 1; }
  sleep 120
done
OMP_NUM_THREADS=4 Rscript logs/dissection/dissection_report.R > logs/dissection/dissection_report.log 2>&1 &&
  echo "$(date '+%F %T') dissection report done (data/02_annotated_data/mouse_dissection_composition.csv, figures/colon/spatial/mouse_dissection.png)" >> $STATUS ||
  echo "$(date '+%F %T') FAILED dissection report (logs/dissection/dissection_report.log)" >> $STATUS
