#!/usr/bin/env bash
# EBCC section 3b-5 (../EBCC_PIPELINE.md, frac by marker AUROC) for colon, alongside 03/04:
#   a. full EBCC fits (M fixed, R = 10) for frac = 0.25 .. 5 (frac = 0 reused from the final fit)
#   b. meanwhile: wait for the 03 run already in progress (pid $1), then all of 04
#   c. choose frac per species; if changed: final fit from the saved fit, 02 export, 03 again
#      (fingerprinted), then only scANVI if cells/genes are unchanged, otherwise all of 04
#   d. all of 05, then figures (colon only)
# Run it in screen:  screen -dmS colon_frac bash run_colon_frac.sh <pid of the running 03 R process>
set -u
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
G=logs/ebcc_frac
STATUS=logs/pipeline_status.txt
FIT_CORES=${FIT_CORES:-4}
stage() { echo "$(date '+%F %T') $*" | tee -a $STATUS; }
trap 'stage "frac runner killed by a signal"; exit 1' TERM INT HUP
R03=$1

stage "EBCC frac fits started (M fixed, R = 10, frac 0.25-5; frac 0 reused)"
(
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 N_THREADS=1   # restarts run as processes
  for s in human mouse; do for f in 0.25 0.5 1 2 5; do echo "$s $f"; done; done |
    xargs -P 10 -L 1 bash -c 'nice -n 5 Rscript '$G'/ebcc_frac_fit.R "$0" "$1" '$FIT_CORES' > '$G'/fit_$0_$1.log 2>&1 || echo "FAILED $0 frac=$1"' \
    > $G/failed.txt
) & FITS=$!

(
  while kill -0 $R03 2>/dev/null; do sleep 60; done
  # the 03 process belonged to the stopped runner, so check its outputs instead of its exit status
  if grep -q '^Error' logs/03_concatenate.log || [ "$(grep -c 'Dropped' logs/03_concatenate.log)" != 4 ] ||
     [ ! -s data/genes8/03_concatenated_data.h5ad ]; then
    stage "FAILED 03 (logs/03_concatenate.log)"; exit 1
  fi
  stage "03 done"
  N_JOBS=4 N_THREADS=6 bash run_04.sh > logs/run_04.log 2>&1
  stage "04 done: $(grep -vc 'exit=0' logs/04/status.txt) failed jobs (logs/04/status.txt)"
) & P04=$!

wait $FITS
fits_ok=1; [ -s $G/failed.txt ] && { fits_ok=0; stage "FAILED EBCC frac fits: $(tr '\n' ' ' < $G/failed.txt)(logs in $G/)"; }
[ $fits_ok = 1 ] && stage "EBCC frac fits done"
wait $P04 || { stage "stopping: 03 failed"; exit 1; }
[ $fits_ok = 1 ] || exit 1

changed=()
for s in human mouse; do
  Rscript $G/ebcc_frac_choose.R $s > $G/choose_$s.log 2>&1 || { stage "FAILED choosing $s frac ($G/choose_$s.log)"; exit 1; }
  stage "$(grep 'chosen frac' $G/choose_$s.log)"
  grep -q '^CHANGED yes' $G/choose_$s.log && changed+=($s)
done

if [ ${#changed[@]} -gt 0 ]; then
  stage "final EBCC fit from the saved fits: ${changed[*]}"
  fail=0
  for s in "${changed[@]}"; do Rscript $G/ebcc_frac_finish.R $s > $G/finish_$s.log 2>&1 & done
  for j in $(jobs -p); do wait $j || fail=1; done
  [ $fail = 1 ] && { stage "FAILED final EBCC fit ($G/finish_*.log)"; exit 1; }
  Rscript export_02_h5ad.R > logs/export_02_h5ad.log 2>&1 || { stage "FAILED 02 export/check (logs/export_02_h5ad.log)"; exit 1; }
  stage "final EBCC fits saved, 02 h5ad exported"

  ${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python $G/fingerprint_03.py > $G/fingerprint_before.txt 2>&1 ||
    { stage "FAILED fingerprint ($G/fingerprint_before.txt)"; exit 1; }
  Rscript -e 'knitr::purl("03_concatenate.qmd", output = "logs/03_concatenate.R", quiet = TRUE, documentation = 0)' &&
    N_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 Rscript logs/03_concatenate.R > logs/03_concatenate_rerun.log 2>&1 ||
    { stage "FAILED 03 rerun (logs/03_concatenate_rerun.log)"; exit 1; }
  ${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python $G/fingerprint_03.py > $G/fingerprint_after.txt 2>&1 ||
    { stage "FAILED fingerprint ($G/fingerprint_after.txt)"; exit 1; }

  if cmp -s $G/fingerprint_before.txt $G/fingerprint_after.txt; then
    stage "03 rerun done, cells and genes unchanged: rerunning scANVI only"
    export N_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
    scanvi() {   # one GPU lane
      local gpu=$1; shift
      for d in "$@"; do
        CUDA_VISIBLE_DEVICES=$gpu ${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python 04.08_integrate_scANVI.py $d \
          > logs/04/04.08_integrate_scANVI__$d.log 2>&1
        echo "04.08_integrate_scANVI.py $d exit=$?" >> logs/04/status_scanvi_rerun.txt
      done
    }
    : > logs/04/status_scanvi_rerun.txt
    scanvi 0 full_mouse genes2 genes6 & scanvi 1 dissected_mouse genes4 genes8 & wait
    stage "scANVI rerun done: $(grep -vc 'exit=0' logs/04/status_scanvi_rerun.txt) failed (logs/04/status_scanvi_rerun.txt)"
  else
    stage "03 rerun done, cells or genes CHANGED ($G/fingerprint_*.txt): rerunning all of 04"
    N_JOBS=4 N_THREADS=8 bash run_04.sh > logs/run_04.log 2>&1
    stage "04 done: $(grep -vc 'exit=0' logs/04/status.txt) failed jobs (logs/04/status.txt)"
  fi
else
  stage "frac unchanged for both species: labels as fitted, no rerun needed"
fi

mkdir -p logs/05; : > logs/05/status.txt
N_JOBS=6 N_THREADS=4 bash run_05.sh > logs/run_05.log 2>&1
stage "05 done: $(grep -vc 'exit=0' logs/05/status.txt) failed jobs (logs/05/status.txt)"

# figures (shared ../analysis/visualisations.R, colon only): fresh UMAP and spatial folders, since it skips existing files
rm -rf ../figures/colon/umap ../figures/colon/spatial
(cd .. && OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 N_THREADS=8 Rscript analysis/visualisations.R colon) > logs/visualisations.log 2>&1 ||
  { stage "FAILED figures (logs/visualisations.log)"; exit 1; }
stage "figures done (figures/colon/)"
