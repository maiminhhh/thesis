#!/usr/bin/env bash
# Label sensitivity check (../EBCC_PIPELINE.md 9b), steps 1-2 only: after the main pipeline (05 with EBCC labels +
# figures, which write marker_auroc.csv) has finished, pick the comparison labels per species. The 05 rerun is
# started separately, after the choice has been reported.
set -u
cd "$(dirname "$0")/.."
until grep -qE "PIPELINE FINISHED|STOP" logs/pipeline.log 2>/dev/null; do sleep 300; done
if grep -q STOP logs/pipeline.log; then echo "STOP: the main pipeline stopped ($(date))"; exit 1; fi
A=../figures/breast_cancer/annotation/marker_auroc.csv
if [ ! -f $A ] || [ $A -ot logs/03/partA.log ]; then echo "STOP: no fresh $A (did the figures step fail?) ($(date))"; exit 1; fi
nice -n 10 Rscript 05_label_sensitivity_choice.R > logs/05_label_sensitivity_choice.log 2>&1 \
  || { echo "STOP: choice failed, see logs/05_label_sensitivity_choice.log"; exit 1; }
tail -4 logs/05_label_sensitivity_choice.log
echo "== CHOICE READY ($(date))"
