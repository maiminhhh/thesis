#!/usr/bin/env bash
# Finishing steps after the early-results reshuffle of 2026-10-26 (logs/finish.log):
#  1. full_mouse main scores (bc_05full) + scVI-latent scores of the other datasets (bc_05scvi) done -> final figures
#  2. (sensitivity rerun on full_mouse: started separately, screen bc_sens_full)
#  3. all sensitivity scores done (bc_sens + full_mouse) -> 05_label_sensitivity.R (tau, top-3, figure)
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/breast_cancer
L=logs/finish.log; log() { echo "$(date '+%F %T') $*" >> $L; }
done_ok() { grep -q '^exit' $1 2>/dev/null; }
LC="hsapiens=celltype_seurat,mmusculus=celltype_seurat_cross"
log "== finish.sh started"
until done_ok logs/05/run_05full.log; do sleep 120; done
log "full_mouse main scores: $(grep -E '^exit' logs/05/run_05full.log) $(grep FAILED logs/05/run_05full.log)"
# (the full_mouse sensitivity rerun was started separately at 19:40 in screen bc_sens_full, alongside the main scores)
until done_ok logs/05/run_05scvi.log; do sleep 120; done
log "scVI-latent scores: $(grep -E '^exit' logs/05/run_05scvi.log)"
( cd .. && OMP_NUM_THREADS=4 nice -n 10 Rscript analysis/visualisations.R breast_cancer ) > logs/figs_final.log 2>&1
log "== FIGURES DONE (exit $?)"
wait
until done_ok logs/05/run_05_sens.log && done_ok logs/05/run_05_sens_full.log; do sleep 120; done
log "sensitivity scores: $(grep -hE '^exit' logs/05/run_05_sens.log logs/05/run_05_sens_full.log | tr '\n' ' ')"
nice -n 10 Rscript 05_label_sensitivity.R > logs/05_label_sensitivity.log 2>&1
log "== SENSITIVITY DONE (exit $?)"
bash pipeline/ntfy_send.sh "breast_cancer" "All done: final figures and label sensitivity (see logs/finish.log)" 3 > /dev/null 2>&1
