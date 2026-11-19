#!/usr/bin/env bash
# After the genes6 rerun (PCA-initialised UMAP): 04 -> 05 -> figures on the 03 datasets on disk, then the
# label-sensitivity choice. Do not edit while it runs.
set -u
cd "$(dirname "$0")/.."
until grep -qE "^exit [0-9]+" logs/03/partB_0.6.log 2>/dev/null; do sleep 120; done
if ! grep -q "^exit 0" logs/03/partB_0.6.log; then echo "STOP: genes6 rerun failed, see logs/03/partB_0.6.log ($(date))"; exit 1; fi
echo "genes6 done: $(grep -E 'removed [0-9]+ of' logs/03/partB_0.6.log) ($(date))"
screen -S bc_sensitivity_choice -X quit 2>/dev/null
screen -dmS bc_pipeline bash -c 'cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/breast_cancer && SKIP_WAIT=1 SKIP_03=1 bash run_pipeline.sh > logs/pipeline.log 2>&1'
sleep 30
screen -dmS bc_sensitivity_choice bash -c 'cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/breast_cancer && bash pipeline/run_sensitivity_choice.sh > logs/sensitivity_choice.log 2>&1'
echo "started bc_pipeline (SKIP_03=1) and bc_sensitivity_choice ($(date))"
