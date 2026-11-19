#!/usr/bin/env bash
# rerun the SAMap jobs that ran out of GPU memory (GPU 1 was shared with fetal_liver's SAMap), one at a time on GPU 1
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/colon
export N_THREADS=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
for d in dissected_mouse genes4 genes8; do
  start=$SECONDS
  CUDA_VISIBLE_DEVICES=1 timeout 12h ${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/samap/bin/python 04.02_integrate_SAMap.py $d > logs/04/04.02_integrate_SAMap__$d.log 2>&1
  rc=$?
  printf '%-40s %-16s exit=%d  %5d min\n' 04.02_integrate_SAMap.py $d $rc $(( (SECONDS - start) / 60 )) >> logs/04/status_samap_rerun.txt
  echo "$(date '+%F %T') SAMap rerun $d: exit=$rc" >> logs/pipeline_status.txt
done
