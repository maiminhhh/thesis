#!/usr/bin/env bash
# 05 rerun after the PCR_batch baseline fix (../EBCC_PIPELINE.md section 7, joint_pca):
#   1. all of 05 (EBCC labels)             2. this tissue's figures (alongside step 3)
#   3. the 9b best-annotator 05 rerun, then label_sensitivity.R
# Run in screen:  screen -dmS fetal_pcr bash rerun_05_pcr_fix.sh
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/fetal_liver
LOG=logs/rerun_downstream.log
say() { echo "$(date '+%F %T') [pcr] $*" | tee -a $LOG; }
fail() { say "STOPPED: $*"; exit 1; }
N=${N_JOBS:-9}

rm -rf logs/05
say "05 started (joint_pca baseline), $N jobs at a time"
N_JOBS=$N N_THREADS=4 bash run_05.sh > logs/run_05.log 2>&1
failed=$(grep -c 'exit=[1-9]' logs/05/status.txt)
[ "$failed" -gt 0 ] && fail "$failed 05 job(s) failed (logs/05/status.txt)"
say "05 done ($(grep -c 'exit=0' logs/05/status.txt) jobs)"

rm -rf ../figures/fetal_liver/umap ../figures/fetal_liver/spatial
say "figures started"
( (cd .. && Rscript analysis/visualisations.R fetal_liver > fetal_liver/logs/visualisations.log 2>&1) \
    && say "figures done: ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/figures/fetal_liver/" \
    || say "STOPPED: visualisations.R failed (logs/visualisations.log)" ) & figs=$!

LABEL_COLS=$(python3 -c "import json; print(json.load(open('data/05_evaluated_score/label_sensitivity_choice.json'))['label_cols'])")
rm -rf logs/05_best_annotator
say "best-annotator 05 started (LABEL_COLS=$LABEL_COLS)"
LABEL_SET=best_annotator LABEL_COLS=$LABEL_COLS N_JOBS=$N N_THREADS=4 bash run_05.sh > logs/run_05_best_annotator.log 2>&1
failed=$(grep -c 'exit=[1-9]' logs/05_best_annotator/status.txt)
[ "$failed" -gt 0 ] && fail "$failed best-annotator 05 job(s) failed (logs/05_best_annotator/status.txt)"
say "best-annotator 05 done ($(grep -c 'exit=0' logs/05_best_annotator/status.txt) jobs)"
Rscript label_sensitivity.R > logs/label_sensitivity_R.log 2>&1 || fail "label_sensitivity.R failed (logs/label_sensitivity_R.log)"
wait $figs
say "all done: scores, figures and label sensitivity"
