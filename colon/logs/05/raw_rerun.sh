#!/usr/bin/env bash
# rescore Unintegrated (05.00_evaluate_raw.py) on every dataset: the 12x2 05 restart stopped its own copies and the
# originals ended with the old runner; one job per dataset, 2 threads each, alongside the main 05 run
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
export N_THREADS=2 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 NUMBA_NUM_THREADS=2
: > logs/05/status_raw_rerun.txt
for d in full_mouse dissected_mouse genes2 genes4 genes6 genes8; do
  ( start=$SECONDS
    ${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/scib/bin/python 05.00_evaluate_raw.py $d > logs/05/05.00_evaluate_raw__${d}_rerun.log 2>&1
    printf '%-40s %-16s exit=%d  %5d min\n' 05.00_evaluate_raw.py $d $? $(( (SECONDS - start) / 60 )) >> logs/05/status_raw_rerun.txt ) &
done
wait
echo "$(date '+%F %T') Unintegrated 05 rescoring done: $(grep -vc 'exit=0' logs/05/status_raw_rerun.txt) failed (logs/05/status_raw_rerun.txt)" >> logs/pipeline_status.txt
