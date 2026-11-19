#!/usr/bin/env bash
# Rescore the graph (knn) outputs that were scored before the utils.score() graph fix (2026-10-25 ~11:10):
# dissected_mouse SAMap + UINMF-BBKNN, and genes2 SAMap once its pre-fix run has exited.
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/breast_cancer
PY=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/scib/bin/python
export N_CPU=8 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
L=logs/05/rescore_knn.log
run() { nice -n 10 $PY $1 $2 > logs/05/${1%.py}_$2.rescore.log 2>&1; echo "$(date '+%F %T') $1 $2: exit $?" >> $L; }
echo "$(date '+%F %T') start" >> $L
run 05.02_evaluate_SAMap.py dissected_mouse &
run 05.17_evaluate_uinmf_BBKNN.py dissected_mouse &
while pgrep -u $USER -f '05.02_evaluate_SAMap.py genes2' | xargs -r -I{} readlink /proc/{}/cwd | grep -qx "$PWD"; do sleep 60; done
run 05.02_evaluate_SAMap.py genes2 &
wait
echo "$(date '+%F %T') all done" >> $L
