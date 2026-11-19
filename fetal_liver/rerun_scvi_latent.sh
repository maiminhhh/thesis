#!/usr/bin/env bash
# scVI-latent methods (../EBCC_PIPELINE.md 6b), alongside rerun_05_pcr_fix.sh. 04.18-04.24 are done;
# the best-annotator scores of all methods (incl. the 7 new ones) run in screen fetal_best
# (JOB_LIST=logs/05_best_annotator_jobs.txt). This chain:
#   1. after the main 05: the 7 new 05 scores with the EBCC labels
#   2. after the best-annotator run and rerun_05_pcr_fix.sh: label_sensitivity.R and the
#      figures again, with all 26 methods
# Run in screen:  screen -dmS fetal_scvi_chain bash rerun_scvi_latent.sh
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/fetal_liver
LOG=logs/rerun_downstream.log
say() { echo "$(date '+%F %T') [scvi] $*" | tee -a $LOG; }
fail() { say "STOPPED: $*"; exit 1; }

# the main 05 (rerun_05_pcr_fix.sh) finished on 2026-10-26 14:03 (114 jobs, 0 failed)
grep -q '\[pcr\] 05 done' $LOG || fail "the main 05 has not finished"
say "05 scVI-latent (EBCC labels) started"
SCRIPTS='05.1[89]*.py 05.2*.py' N_JOBS=9 N_THREADS=4 bash run_05.sh > logs/run_05_scvi.log 2>&1
n=$(grep -E '^05\.(1[89]|2[0-4])_' logs/05/status.txt | grep -c 'exit=[1-9]')
[ "$n" -gt 0 ] && fail "$n scVI-latent 05 job(s) failed (logs/05/status.txt)"
say "05 scVI-latent (EBCC labels) done"

say "waiting for the best-annotator 05 and rerun_05_pcr_fix.sh"
while kill -0 "$(head -1 logs/05_best_annotator.lock)" 2>/dev/null || screen -ls | grep -q '\.fetal_pcr'; do sleep 300; done
n=$(grep -c 'exit=[1-9]' logs/05_best_annotator_early/status.txt)
[ "$n" -gt 0 ] && fail "$n best-annotator 05 job(s) failed (logs/05_best_annotator_early/status.txt)"
Rscript label_sensitivity.R > logs/label_sensitivity_R.log 2>&1 || fail "label_sensitivity.R failed (logs/label_sensitivity_R.log)"
say "label sensitivity redone with all 26 methods"
(cd .. && Rscript analysis/visualisations.R fetal_liver > fetal_liver/logs/visualisations.log 2>&1) || fail "visualisations.R failed (logs/visualisations.log)"
# (the original run copied the separately made mouse_dissection.png back into ../figures/fetal_liver/spatial/ here)
say "all done: 26 methods scored, figures and label sensitivity redone"
