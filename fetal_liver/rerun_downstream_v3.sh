#!/usr/bin/env bash
# Fallback for rerun_downstream_v2.sh when its only failure is SAMap full_mouse hitting the
# 8h job timeout (exit 124). The EBCC labels did not change (frac selection: human M = 5,
# frac = 0; mouse M = 3, frac = 0), so 03 is not rerun.
#   1. wait for v2 to exit; do nothing unless it stopped on SAMap full_mouse alone
#   2. SAMap full_mouse again (GPU1, 20h limit), alongside:
#   3. scANVI on every dataset, then 05 on everything except SAMap full_mouse
#   4. when SAMap is done, 05.02 on full_mouse; then this tissue's figures
# Run in screen:  screen -dmS fetal_rerun3 bash rerun_downstream_v3.sh
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/fetal_liver
LOG=logs/rerun_downstream.log
CONDA=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}
export N_THREADS=4 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 \
       NUMEXPR_NUM_THREADS=4 NUMBA_NUM_THREADS=4
say() { echo "$(date '+%F %T') [v3] $*" | tee -a $LOG; }
fail() { say "STOPPED: $*"; exit 1; }

# 1. ---- v2 ----
while pgrep -f 'bash rerun_downstream_v2.sh' > /dev/null; do sleep 60; done
tail -1 $LOG | grep -q '\[v2\] STOPPED: 1 04 job(s) failed' || { echo "v2 did not stop on a single 04 failure; v3 not needed"; exit 0; }
failed=$(tail -n "+$(grep -n '^# run ' logs/04/status.txt | tail -1 | cut -d: -f1)" logs/04/status.txt | grep 'exit=[1-9]')
echo "$failed" | grep -q '^04.02_integrate_SAMap.py *full_mouse *exit=124' && [ $(echo "$failed" | wc -l) -eq 1 ] \
  || fail "v2's failed 04 job is not the SAMap full_mouse timeout: $failed"
# every other 04 job of the rerun succeeded (all '# run' sections since v2 started)
others=$(awk '/^# run 2026-10-23 14:46/ {on = 1} on' logs/04/status.txt | grep 'exit=[1-9]' | grep -v '^04.02_integrate_SAMap.py *full_mouse *exit=124')
[ -z "$others" ] || fail "other failed 04 jobs: $others"
say "v2 stopped on the SAMap full_mouse timeout only; taking over"

# 2. ---- SAMap full_mouse, 20h limit ----
samap_job() {
  local start=$SECONDS
  CUDA_VISIBLE_DEVICES=1 timeout 20h $CONDA/samap/bin/python -u 04.02_integrate_SAMap.py full_mouse \
    > logs/04/04.02_integrate_SAMap__full_mouse.log 2>&1
  local rc=$?
  printf '%-40s %-16s exit=%d  %5d min\n' 04.02_integrate_SAMap.py full_mouse $rc $(( (SECONDS - start) / 60 )) >> logs/04/status.txt
  return $rc
}
say "SAMap full_mouse restarted (20h limit)"
( samap_job ) & samap=$!

# 3. ---- scANVI, then 05 without SAMap full_mouse ----
say "scANVI started"
N_JOBS=2 ONLY='04.08_' bash run_04.sh > logs/run_04_scanvi.log 2>&1
failed=$(tail -n "+$(grep -n '^# run ' logs/04/status.txt | tail -1 | cut -d: -f1)" logs/04/status.txt | grep -c 'exit=[1-9]')
[ "$failed" -gt 0 ] && fail "$failed scANVI job(s) failed (logs/04/status.txt)"
say "scANVI done"
rm -rf logs/05
say "05 started (SAMap full_mouse scored when its 04 rerun finishes)"
N_JOBS=6 N_THREADS=4 bash run_05.sh dissected_mouse genes2 genes4 genes6 genes8 > logs/run_05.log 2>&1
N_JOBS=6 N_THREADS=4 SCRIPTS="$(ls 05.*.py | grep -v '^05.02_')" bash run_05.sh full_mouse >> logs/run_05.log 2>&1

# 4. ---- SAMap full_mouse: 05, figures ----
say "waiting for SAMap full_mouse"
wait $samap || fail "SAMap full_mouse failed again (logs/04/04.02_integrate_SAMap__full_mouse.log)"
N_JOBS=1 N_THREADS=4 SCRIPTS=05.02_evaluate_SAMap.py bash run_05.sh full_mouse >> logs/run_05.log 2>&1
failed=$(grep -c 'exit=[1-9]' logs/05/status.txt)
[ "$failed" -gt 0 ] && fail "$failed 05 job(s) failed (logs/05/status.txt)"
say "05 done ($(grep -c 'exit=0' logs/05/status.txt) jobs)"
rm -rf ../figures/fetal_liver/umap ../figures/fetal_liver/spatial
say "figures started"
(cd .. && Rscript analysis/visualisations.R fetal_liver > fetal_liver/logs/visualisations.log 2>&1) \
  || fail "visualisations.R failed (logs/visualisations.log)"
say "all done: figures in ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/figures/fetal_liver/"
