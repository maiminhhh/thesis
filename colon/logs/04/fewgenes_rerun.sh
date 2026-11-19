#!/usr/bin/env bash
# rerun the jobs that failed for too few shared genes, after capping LIGER k / RPCA-UINMF dims
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
export N_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
for job in "04.10_integrate_LIGER.R genes6" "04.10_integrate_LIGER.R genes8" "04.16_integrate_uinmf_RPCA.R genes8"; do
  set -- $job; start=$SECONDS
  timeout 6h Rscript $1 $2 > logs/04/${1%.*}__$2.log 2>&1; rc=$?
  printf '%-40s %-16s exit=%d  %5d min\n' $1 $2 $rc $(( (SECONDS - start) / 60 )) >> logs/04/status_fewgenes_rerun.txt
  echo "$(date '+%F %T') few-genes rerun $1 $2: exit=$rc" >> logs/pipeline_status.txt
done
