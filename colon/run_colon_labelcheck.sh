#!/usr/bin/env bash
# Targeted label sensitivity check (section 9b of ../EBCC_PIPELINE.md), in two steps:
#   choose:  wait for run_colon_frac.sh's "figures done" (marker_auroc.csv), then label_choice.R
#            -> data/05_evaluated_score/label_sensitivity_choice.json; stops there (the choice is reported first)
#   rescore: rescore all of 05 with the best-annotator labels, then label_sensitivity.R
# Run in screen:  screen -dmS colon_labels bash run_colon_labelcheck.sh choose
#                 screen -dmS colon_labels bash run_colon_labelcheck.sh rescore
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
STATUS=logs/pipeline_status.txt
CHOICE=data/05_evaluated_score/label_sensitivity_choice.json
stage() { echo "$(date '+%F %T') $*" | tee -a $STATUS; }
trap 'stage "label check runner ($1) killed by a signal"; exit 1' TERM INT HUP

case $1 in
choose)
  n=$(wc -l < $STATUS)
  until tail -n +$((n + 1)) $STATUS | grep -q 'figures done'; do
    screen -ls | grep -q colon_frac || { stage "label check: colon_frac exited without drawing the figures, stopping"; exit 1; }
    sleep 120
  done
  Rscript label_choice.R > logs/label_choice.log 2>&1 || { stage "FAILED label choice (logs/label_choice.log)"; exit 1; }
  stage "label check choice: $(grep -E '^(human|mouse):' logs/label_choice.log | tr '\n' ' ')-> $(tail -1 logs/label_choice.log)"
  ;;
rescore)
  cols=$(sed -n 's/^LABEL_COLS=//p' logs/label_choice.log)
  [ -n "$cols" ] || { stage "label check: no LABEL_COLS in logs/label_choice.log (EBCC best in both species?), nothing to rescore"; exit 1; }
  stage "05 rescoring with best-annotator labels started ($cols)"
  rm -rf logs/05_best_annotator
  LABEL_SET=best_annotator LABEL_COLS=$cols N_JOBS=4 N_THREADS=6 bash run_05.sh > logs/run_05_best_annotator.log 2>&1
  stage "05 best-annotator rescoring done: $(grep -vc 'exit=0' logs/05_best_annotator/status.txt) failed jobs (logs/05_best_annotator/status.txt)"
  Rscript label_sensitivity.R > logs/label_sensitivity.log 2>&1 || { stage "FAILED label sensitivity (logs/label_sensitivity.log)"; exit 1; }
  stage "label sensitivity done (data/05_evaluated_score/label_sensitivity*.csv, figures/colon/scores/label_sensitivity.png)"
  ;;
*) echo "usage: bash run_colon_labelcheck.sh choose|rescore"; exit 1 ;;
esac
