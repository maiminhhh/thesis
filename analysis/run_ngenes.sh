#!/bin/bash
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}
LOG=${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/figures/orthologue_ngenes_run.log
echo "===== analysis/orthologue_input.R (genes model with BBKNN)" > $LOG
if R_LIBS_SITE= R_LIBS_USER= ${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/r45/bin/Rscript --vanilla analysis/orthologue_input.R >> $LOG 2>&1; then
  echo "===== ALL DONE" >> $LOG
else
  echo "===== FAILED: analysis/orthologue_input.R" >> $LOG
fi
