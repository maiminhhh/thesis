#!/usr/bin/env bash
# Smoke test of the scVI-latent methods (../EBCC_PIPELINE.md 6b) on one dataset: 04.18 on GPU 0, then 04.19-04.24,
# then the 05.18-05.24 scores from pipeline/staging_05scvi/ (kept out of the project folder while the main 05 runs,
# whose dispatcher globs 05*.py). Log lines "exit N" go to logs/04scvi/smoke.log.
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/breast_cancer
g=${1:-genes8}; T=${THREADS:-4}
export RETICULATE_PYTHON=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python
export N_THREADS=$T N_CPU=$T OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T NUMEXPR_NUM_THREADS=$T
PY_TORCH=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/torch/bin/python
PY_SCIB=${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/scib/bin/python
L=logs/04scvi; mkdir -p $L; S=$L/smoke.log
job() { local name=$1; shift; nice -n 10 "$@" > $L/${name}_$g.log 2>&1; echo "$(date '+%F %T') $name $g: exit $?" >> $S; }
echo "$(date '+%F %T') smoke start $g" >> $S
CUDA_VISIBLE_DEVICES=0 job 04.18 $PY_TORCH 04.18_scVI_latent.py $g
[ -s data/$g/scvi_latent.csv ] || { echo "$(date '+%F %T') STOP: no latent" >> $S; exit 1; }
for f in 04.19_integrate_fastMNN_scvi 04.20_integrate_CCA_scvi 04.21_integrate_Harmony_scvi 04.23_integrate_RPCA_scvi; do
  CUDA_VISIBLE_DEVICES="" job ${f%%_*} Rscript $f.R $g & done
CUDA_VISIBLE_DEVICES="" job 04.22 $PY_TORCH 04.22_integrate_Scanorama_scvi.py $g &
CUDA_VISIBLE_DEVICES="" job 04.24 $PY_TORCH 04.24_integrate_BBKNN_scvi.py $g &
wait
for f in pipeline/staging_05scvi/05.*.py; do
  b=$(basename $f .py); PYTHONPATH=. job ${b%%_*} $PY_SCIB $f $g &
  while [ $(jobs -r | wc -l) -ge 4 ]; do sleep 30; done
done
wait
echo "$(date '+%F %T') smoke all done $g" >> $S
