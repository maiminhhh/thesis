#!/usr/bin/env bash
# After the main 05 (rerun_downstream_v2.sh): rescore 05 once with the best-annotator labels
# (../EBCC_PIPELINE.md 9b; LABEL_COLS from data/05_evaluated_score/label_sensitivity_choice.json),
# then label_sensitivity.R.
# Run in screen:  screen -dmS fetal_labelsens bash label_sensitivity_chain.sh
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/fetal_liver
LOG=logs/label_sensitivity.log
say() { echo "$(date '+%F %T') $*" | tee -a $LOG; }
LABEL_COLS=$(python3 -c "import json; print(json.load(open('data/05_evaluated_score/label_sensitivity_choice.json'))['label_cols'])")
say "waiting for the main 05 (logs/rerun_downstream.log); then LABEL_COLS=$LABEL_COLS"
# the main 05 is run by rerun_downstream_v2.sh, or by v3 if v2 stopped on the SAMap full_mouse timeout
until grep -q '\[v[23]\] 05 done' logs/rerun_downstream.log; do
  grep -q '\[v3\] STOPPED' logs/rerun_downstream.log && { say "STOPPED: the v3 driver stopped before 05 finished"; exit 1; }
  grep -q '\[v2\] STOPPED' logs/rerun_downstream.log && ! pgrep -f 'bash rerun_downstream_v3.sh' > /dev/null \
    && { say "STOPPED: the v2 driver stopped before 05 finished (no v3 running)"; exit 1; }
  sleep 120
done
say "05 with best-annotator labels started"
rm -rf logs/05_best_annotator
LABEL_SET=best_annotator LABEL_COLS=$LABEL_COLS N_JOBS=6 N_THREADS=4 bash run_05.sh > logs/run_05_best_annotator.log 2>&1
failed=$(grep -c 'exit=[1-9]' logs/05_best_annotator/status.txt)
[ "$failed" -gt 0 ] && { say "STOPPED: $failed best-annotator 05 job(s) failed (logs/05_best_annotator/status.txt)"; exit 1; }
say "05 with best-annotator labels done ($(grep -c 'exit=0' logs/05_best_annotator/status.txt) jobs)"
Rscript label_sensitivity.R > logs/label_sensitivity_R.log 2>&1 || { say "STOPPED: label_sensitivity.R failed (logs/label_sensitivity_R.log)"; exit 1; }
say "all done: data/05_evaluated_score/label_sensitivity*.csv, ../figures/fetal_liver/scores/label_sensitivity.png"
