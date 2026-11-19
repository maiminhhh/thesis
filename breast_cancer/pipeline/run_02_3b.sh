#!/usr/bin/env bash
# Section 3b of ../EBCC_PIPELINE.md for one species: M from the held-out grid (1-SE rule on the frac = 0
# column), six full EBCC fits (frac = 0, 0.25, 0.5, 1, 2, 5; R = 10), then frac by marker AUROC and the final save.
#   bash pipeline/run_02_3b.sh <human|mouse> [launch]
# With "launch" the six fits are started here (4 threads each); without it, they are expected to be
# running already and the script only waits for them. Run it in screen. Do not edit while it runs.
set -u
cd "$(dirname "$0")/.."
sp=$1; launch=${2:-}
D=data/02_annotated_data
step() { echo "== $* ($(date))"; }

step "$sp: waiting for the held-out grid"
until [ -f $D/${sp}_ebcc_heldout_grid_chosen.csv ] && ! screen -ls | grep -q "bc_ebcc_$sp"; do sleep 120; done
M=$(Rscript -e "source('ebcc_grid.R'); cat(choose_M_frac0(read.csv('$D/${sp}_ebcc_heldout_grid.csv')))")
step "$sp: M = $M (1-SE rule on the frac = 0 column)"

fracs="0 0.25 0.5 1 2 5"
if [ "$launch" = launch ]; then
  export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 N_THREADS=4
  for f in $fracs; do
    nice -n 5 Rscript pipeline/02_${sp}_ebcc_frac_fit.R $M $f > logs/02/${sp}_frac_$f.log 2>&1 &
  done
  step "$sp: launched the six frac fits"
fi
for f in $fracs; do
  until [ -f $D/${sp}_ebcc_frac_fits/M${M}_frac$f.rds ]; do
    if grep -qE "Error|Execution halted" logs/02/${sp}_frac_$f.log 2>/dev/null; then
      echo "STOP: $sp fit frac = $f failed, see logs/02/${sp}_frac_$f.log"; exit 1
    fi
    sleep 120
  done
  step "$sp: fit frac = $f saved"
done
wait

step "$sp: frac by marker AUROC + final save"
nice -n 5 Rscript pipeline/02_${sp}_ebcc_frac_select.R > logs/02/${sp}_frac_select.log 2>&1 \
  || { echo "STOP: $sp frac selection failed, see logs/02/${sp}_frac_select.log"; exit 1; }
grep -E "^$sp: M = |final fit|unchanged" logs/02/${sp}_frac_select.log
cat $D/${sp}_ebcc_frac_auroc.csv
step "$sp: 3b FINISHED"
