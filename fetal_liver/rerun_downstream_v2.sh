#!/usr/bin/env bash
# Downstream rerun, second version (../EBCC_PIPELINE.md as of 2026-10-23 13:39: frac chosen
# by marker AUROC). Replaces rerun_downstream.sh, which stops itself after its 03 run.
#   1. wait for that 03 run; check cells / genes against the fingerprint taken before it
#   2. label-independent 04, alongside the frac fits: eval_cells.txt (full / dissected);
#      genes2-8 everything except SAMap (raw counts only) and scANVI (labels);
#      full / dissected scVI (saves the model scANVI starts from); SAMap full_mouse
#   3. wait for the frac selection; if labels changed, rerun 03 and check cells / genes again
#   4. scANVI on every dataset, then all of 05, then this tissue's figures
# Run in screen:  screen -dmS fetal_rerun2 bash rerun_downstream_v2.sh
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/fetal_liver
LOG=logs/rerun_downstream.log
CONDA=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}
export N_THREADS=4 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 \
       NUMEXPR_NUM_THREADS=4 NUMBA_NUM_THREADS=4
say() { echo "$(date '+%F %T') [v2] $*" | tee -a $LOG; }
fail() { say "STOPPED: $*"; exit 1; }

fingerprint() {
  $CONDA/torch/bin/python -I -c "
import anndata as ad, hashlib, sys
for f in sys.argv[1:]:
    a = ad.read_h5ad(f, backed='r')
    h = lambda x: hashlib.md5('\n'.join(map(str, x)).encode()).hexdigest()[:12]
    print(f, a.n_obs, a.n_vars, h(a.obs_names), h(a.var_names))
" data/03_concatenated_data/03_concatenated_data_{full,subset}.h5ad data/genes{2,4,6,8}/03_concatenated_data.h5ad
}
same_cells() {   # every dataset's cells and genes identical in the two fingerprint files
  diff <(cut -d' ' -f1-5 "$1") <(cut -d' ' -f1-5 "$2") > /dev/null
}

# 1. ---- the 03 run started by rerun_downstream.sh ----
say "waiting for the running 03"
# only this project's 03 (other tissues' pipelines run their own 03_concatenate.R)
fetal_03_running() {
  for p in $(pgrep -f '03_concatenate.R'); do
    [ "$(readlink /proc/$p/cwd 2>/dev/null)" == "$PWD" ] && return 0
  done
  return 1
}
while fetal_03_running; do sleep 60; done
grep -q 'Dropped 25 of 123' logs/03_concatenate.log && ! grep -q '^Error\|Execution halted' logs/03_concatenate.log \
  || fail "03 did not finish cleanly (logs/03_concatenate.log)"
fingerprint > logs/fingerprint_03_after.txt
same_cells logs/fingerprint_03_before_true.txt logs/fingerprint_03_after.txt \
  || fail "cells/genes changed in 03 (logs/fingerprint_03_*): all of 04 needs rerunning; not done automatically"
say "03 done; cells and genes unchanged in all 6 datasets"

# 2. ---- label-independent 04 ----
for d in full_mouse dissected_mouse; do
  $CONDA/torch/bin/python -I -c "
import anndata as ad, pandas as pd, sys
out, concat = sys.argv[1], sys.argv[2]
H = pd.read_csv(out + 'uinmf_H_raw.csv', index_col=0)
ok = set(H.index[(H != 0).sum(axis=1) >= 2])
cells = [c for c in ad.read_h5ad(concat, backed='r').obs_names if c in ok]
open(out + 'eval_cells.txt', 'w').write('\n'.join(cells) + '\n')
print(len(cells))
" data/04_integrated_data/$d/ data/03_concatenated_data/03_concatenated_data_$([ $d == full_mouse ] && echo full || echo subset).h5ad \
    > logs/eval_cells_$d.txt || fail "eval_cells.txt for $d"
  say "eval_cells.txt for $d: $(cat logs/eval_cells_$d.txt) cells"
done
gpu_job() {   # gpu env script dataset
  local gpu=$1 env=$2 s=$3 d=$4 start=$SECONDS
  CUDA_VISIBLE_DEVICES=$gpu timeout 8h $CONDA/$env/bin/python -u $s $d > logs/04/${s%.*}__$d.log 2>&1
  local rc=$?
  printf '%-40s %-16s exit=%d  %5d min\n' $s $d $rc $(( (SECONDS - start) / 60 )) >> logs/04/status.txt
  return $rc
}
say "04 (label-independent) started"
( gpu_job 1 samap 04.02_integrate_SAMap.py full_mouse ) & samap=$!
( gpu_job 0 torch 04.07_integrate_scVI.py full_mouse; gpu_job 0 torch 04.07_integrate_scVI.py dissected_mouse ) & scvi=$!
N_JOBS=4 SKIP='04.0[28]_' bash run_04.sh genes2 genes4 genes6 genes8 > logs/run_04_genes.log 2>&1
wait $samap $scvi
failed=$(tail -n "+$(grep -n '^# run ' logs/04/status.txt | tail -1 | cut -d: -f1)" logs/04/status.txt | grep -c 'exit=[1-9]')
[ "$failed" -gt 0 ] && fail "$failed 04 job(s) failed (logs/04/status.txt)"
say "04 (label-independent) done"

# 3. ---- frac selection; 03 again if labels changed ----
say "waiting for the EBCC frac selection"
until grep -q 'all done' logs/02_ebcc/frac_driver.log; do sleep 120; done
grep -q 'exit=[1-9]' logs/02_ebcc/frac_driver.log && fail "an EBCC frac fit / selection failed (logs/02_ebcc/)"
changed=0
for sp in human mouse; do
  IFS=, read -r _ M frac ch _ _ < <(tail -n 1 data/02_annotated_data/${sp}_ebcc_frac_auroc_chosen.csv)
  sed -i "s/^ebcc_previous = list(M = [0-9.]*, frac = [0-9.]*)/ebcc_previous = list(M = $M, frac = $frac)/" 02_${sp}_annotate.R
  say "$sp: M = $M, frac = $frac by marker AUROC (changed: $ch); recorded in 02_${sp}_annotate.R"
  [ "$ch" == TRUE ] && changed=1
done
if [ $changed == 1 ]; then
  say "labels changed: 03 rerun started"
  Rscript 03_concatenate.R > logs/03_concatenate_final.log 2>&1 || fail "03 failed (logs/03_concatenate_final.log)"
  fingerprint > logs/fingerprint_03_final.txt
  same_cells logs/fingerprint_03_before_true.txt logs/fingerprint_03_final.txt \
    || fail "cells/genes changed in the final 03: all of 04 needs rerunning; not done automatically"
  say "03 rerun done; cells and genes unchanged"
fi

# 4. ---- scANVI, 05, figures ----
say "scANVI started"
N_JOBS=2 ONLY='04.08_' bash run_04.sh > logs/run_04_scanvi.log 2>&1
failed=$(tail -n "+$(grep -n '^# run ' logs/04/status.txt | tail -1 | cut -d: -f1)" logs/04/status.txt | grep -c 'exit=[1-9]')
[ "$failed" -gt 0 ] && fail "$failed scANVI job(s) failed (logs/04/status.txt)"
say "scANVI done"
rm -rf logs/05
say "05 started"
N_JOBS=6 N_THREADS=4 bash run_05.sh > logs/run_05.log 2>&1
failed=$(grep -c 'exit=[1-9]' logs/05/status.txt)
[ "$failed" -gt 0 ] && fail "$failed 05 job(s) failed (logs/05/status.txt)"
say "05 done ($(grep -c 'exit=0' logs/05/status.txt) jobs)"
rm -rf ../figures/fetal_liver/umap ../figures/fetal_liver/spatial
say "figures started"
(cd .. && Rscript analysis/visualisations.R fetal_liver > fetal_liver/logs/visualisations.log 2>&1) \
  || fail "visualisations.R failed (logs/visualisations.log)"
say "all done: figures in ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/figures/fetal_liver/"
