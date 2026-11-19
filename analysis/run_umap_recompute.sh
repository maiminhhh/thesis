#!/bin/bash
cd ${PROJECT_ROOT:?set PROJECT_ROOT to the repository folder}
echo "===== bm_umap_compute $(date)"; ${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/scib/bin/python analysis/bm_umap_compute.py || { echo "===== FAILED: bm_umap_compute"; exit 1; }
echo "===== appendix_unintegrated_umap_compute $(date)"; ${CONDA_ENVS:?set CONDA_ENVS to the folder holding the conda environments}/scib/bin/python analysis/appendix_unintegrated_umap_compute.py || { echo "===== FAILED: unintegrated compute"; exit 1; }
echo "===== ALL DONE $(date)"
