#!/bin/bash
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}
LOG=${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}/figures/mcb_s24_run.log
echo "===== analysis/mcb_analysis.R (section 2.4 model)" > $LOG
if R_LIBS_SITE= R_LIBS_USER= ${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/r45/bin/Rscript --vanilla analysis/mcb_analysis.R >> $LOG 2>&1; then
  echo "===== ALL DONE" >> $LOG
else
  echo "===== FAILED: analysis/mcb_analysis.R" >> $LOG
fi
