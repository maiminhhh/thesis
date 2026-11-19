#!/usr/bin/env bash
# Downstream rerun after the EBCC (M, frac) grid, following ../EBCC_PIPELINE.md, for the
# changes of 2026-10-23: gene-drop datasets renormalised (03), common evaluation cells
# (eval_cells.txt, 05), scANVI starting from 04.07's saved scVI model.
#   1. wait for the EBCC grid (02) and the SAMap queue (they read 02 / 03 outputs)
#   2. write the chosen M / frac back into the 02 scripts
#   3. fingerprint the concatenated data, rerun all of 03, fingerprint again
#   4. 04: genes2-8 everything except SAMap (reads raw counts only, unaffected);
#      full / dissected: scVI + scANVI, SAMap full_mouse, eval_cells.txt from the saved H
#   5. 05: everything      6. figures: visualisations.R fetal_liver
# Run in screen:  screen -dmS fetal_rerun bash rerun_downstream.sh
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/fetal_liver
LOG=logs/rerun_downstream.log
CONDA=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}
export N_THREADS=4 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 \
       NUMEXPR_NUM_THREADS=4 NUMBA_NUM_THREADS=4
say() { echo "$(date '+%F %T') $*" | tee -a $LOG; }
fail() { say "STOPPED: $*"; exit 1; }

# 1. ---- wait for the EBCC grid and the SAMap queue ----
say "waiting for the EBCC grid"
until grep -q 'all done' logs/02_ebcc/driver.log; do sleep 120; done
grep -q 'exit=[1-9]' logs/02_ebcc/driver.log && fail "an EBCC grid / fit job failed (logs/02_ebcc/)"
say "EBCC grid done: $(tail -n +2 data/02_annotated_data/human_ebcc_heldout_grid_chosen.csv) | $(tail -n +2 data/02_annotated_data/mouse_ebcc_heldout_grid_chosen.csv)"
say "waiting for the SAMap queue (genes4/6/8)"
while pgrep -f 'samap_gpu.sh' > /dev/null; do sleep 120; done

# 2. ---- chosen M / frac back into the 02 scripts ----
for sp in human mouse; do
  IFS=, read -r _ M frac changed _ _ < <(tail -n 1 data/02_annotated_data/${sp}_ebcc_heldout_grid_chosen.csv)
  sed -i "s/^ebcc_previous = list(M = [0-9.]*, frac = [0-9.]*)/ebcc_previous = list(M = $M, frac = $frac)/" 02_${sp}_annotate.R
  say "$sp: M = $M, frac = $frac (changed: $changed); recorded in 02_${sp}_annotate.R"
done

# 3. ---- 03, with a fingerprint of cells and genes before and after ----
fingerprint() {
  $CONDA/torch/bin/python -I -c "
import anndata as ad, hashlib, sys
for f in sys.argv[1:]:
    a = ad.read_h5ad(f, backed='r')
    h = lambda x: hashlib.md5('\n'.join(map(str, x)).encode()).hexdigest()[:12]
    print(f, a.n_obs, a.n_vars, h(a.obs_names), h(a.var_names))
" data/03_concatenated_data/03_concatenated_data_{full,subset}.h5ad data/genes{2,4,6,8}/03_concatenated_data.h5ad
}
fingerprint > logs/fingerprint_03_before.txt
say "03 started"
Rscript 03_concatenate.R > logs/03_concatenate.log 2>&1 || fail "03 failed (logs/03_concatenate.log)"
fingerprint > logs/fingerprint_03_after.txt
say "03 done; cells/genes changed in: $(diff <(awk '{print $1,$2,$3,$4,$5}' logs/fingerprint_03_before.txt) \
     <(awk '{print $1,$2,$3,$4,$5}' logs/fingerprint_03_after.txt) | grep '^>' | awk '{print $2}' | tr '\n' ' ')none-other"
# full / dissected reuse their 04 results (only scVI / scANVI rerun), so their cells must not change
for f in data/03_concatenated_data/03_concatenated_data_full.h5ad data/03_concatenated_data/03_concatenated_data_subset.h5ad; do
  [ "$(grep "^$f " logs/fingerprint_03_before.txt)" == "$(grep "^$f " logs/fingerprint_03_after.txt)" ] \
    || fail "cells/genes of $f changed: full/dissected need all of 04 rerun; not done automatically"
done

# 4. ---- 04 ----
# eval_cells.txt for full / dissected from their saved UINMF H (genes2-8 get it from 04.11)
for d in full_mouse dissected_mouse; do
  $CONDA/torch/bin/python -I -c "
import anndata as ad, pandas as pd, sys
out, concat = sys.argv[1], sys.argv[2]
H = pd.read_csv(out + 'uinmf_H_raw.csv', index_col=0)
ok = H.index[(H != 0).sum(axis=1) >= 2]
cells = [c for c in ad.read_h5ad(concat, backed='r').obs_names if c in set(ok)]
open(out + 'eval_cells.txt', 'w').write('\n'.join(cells) + '\n')
print(len(cells))
" data/04_integrated_data/$d/ data/03_concatenated_data/03_concatenated_data_$([ $d == full_mouse ] && echo full || echo subset).h5ad \
    > /tmp/fetal_eval_$d.txt || fail "eval_cells.txt for $d"
  say "eval_cells.txt for $d: $(cat /tmp/fetal_eval_$d.txt) cells"
done
gpu_job() {   # gpu env script dataset
  local gpu=$1 env=$2 s=$3 d=$4 start=$SECONDS
  CUDA_VISIBLE_DEVICES=$gpu timeout 8h $CONDA/$env/bin/python -u $s $d > logs/04/${s%.*}__$d.log 2>&1
  local rc=$?
  printf '%-40s %-16s exit=%d  %5d min\n' $s $d $rc $(( (SECONDS - start) / 60 )) >> logs/04/status.txt
  return $rc
}
say "04 started: genes2-8 (all but SAMap); full/dissected scVI + scANVI; SAMap full_mouse"
( gpu_job 1 samap 04.02_integrate_SAMap.py full_mouse ) & samap=$!
( for d in full_mouse dissected_mouse; do
    gpu_job 0 torch 04.07_integrate_scVI.py $d && gpu_job 0 torch 04.08_integrate_scANVI.py $d
  done ) & scvi=$!
N_JOBS=6 SKIP='04.02_' bash run_04.sh genes2 genes4 genes6 genes8 > logs/run_04_genes.log 2>&1
wait $samap $scvi
failed=$(tail -n "+$(grep -n '^# run ' logs/04/status.txt | tail -1 | cut -d: -f1)" logs/04/status.txt | grep -c 'exit=[1-9]')
[ "$failed" -gt 0 ] && fail "$failed 04 job(s) failed (logs/04/status.txt)"
say "04 done"

# 5. ---- 05 ----
rm -rf logs/05
say "05 started"
N_JOBS=6 N_THREADS=4 bash run_05.sh > logs/run_05.log 2>&1
failed=$(grep -c 'exit=[1-9]' logs/05/status.txt)
[ "$failed" -gt 0 ] && fail "$failed 05 job(s) failed (logs/05/status.txt)"
say "05 done ($(grep -c 'exit=0' logs/05/status.txt) jobs)"

# 6. ---- figures ----
rm -rf ../figures/fetal_liver/umap ../figures/fetal_liver/spatial
say "figures started"
(cd .. && Rscript analysis/visualisations.R fetal_liver > fetal_liver/logs/visualisations.log 2>&1) \
  || fail "visualisations.R failed (logs/visualisations.log)"
say "all done: figures in ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/figures/fetal_liver/"
