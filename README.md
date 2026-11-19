# Cross-species integration benchmark on spatial transcriptomics (human and mouse)

Code, result tables and figures for my thesis. The project benchmarks methods for integrating
human and mouse single-cell spatial transcriptomics data across four tissues (bone marrow, breast cancer,
colon, fetal liver). Cell types are annotated by combining several automatic annotators with an EBCC
consensus model, and each integration is then scored for batch correction and conservation of biology.

> Thesis title, author, supervisor and submission date: _fill in before submission._

## What is in this repository

| Folder / file | Contents |
|---|---|
| `bone_marrow/`, `breast_cancer/`, `colon/`, `fetal_liver/` | One pipeline per tissue, numbered by stage (see below), with the shell scripts that ran it. `data/` inside each holds only the small result tables (scores, model-selection tables, timings). |
| `analysis/` | Cross-tissue statistics, tables and figures for the thesis chapters (index below). |
| `library.R`, `source.R` | R packages and shared functions used by every pipeline (QC, label transfer, EBCC, concatenation on orthologues, integration helpers). |
| `homology_tbl.csv` | Human-mouse orthologue table used to merge the two species. |
| `figures/` | Figures and the tables behind them, per tissue and across tissues. |
| `bootstrap*/`, `sas/`, `extra_metrics/`, `label_map_alt/` | Result tables of the subsampling (bootstrap) analyses and the additional metrics. |
| `docs/` | `EBCC_PIPELINE.md` (the annotation model-selection procedure), `methods_preprocessing_code.md` (pre-processing code excerpts), `thesis_listings/` (code listings shown in the thesis appendix). |
| `envs/` | Conda environment specifications and the R package versions. |

### Pipeline stages (in each tissue folder)

| Stage | Scripts | What it does |
|---|---|---|
| 01 | `01_raw_data_process*` | Load the raw data of each species, quality control |
| 02 | `02_human_annotate*`, `02_mouse_annotate*` | Cell-type annotation: Seurat and SingleR label transfer with same-species and cross-species references, combined by EBCC |
| 03 | `03_concatenate*` | Merge the species on one-to-one orthologues; full and dissected mouse sections; datasets with 20-80% of shared genes removed (`genes2` ... `genes8`) |
| 04 | `04.xx_integrate_<method>` | Integration: FastMNN, SAMap, Seurat CCA and RPCA, Harmony, Scanorama, scVI, scANVI, BBKNN, LIGER, LIGER-UINMF, and the same methods on UINMF and scVI-latent input |
| 05 | `05.xx_evaluate_<method>` | Score each integration (scIB metrics; `utils.py`) |
| 06-08 | bone marrow only for 06/07 | Label transfer, colocalisation, run times |

`run_04.sh`, `run_05.sh` and the other `run_*.sh` / `rerun_*.sh` files are the job scripts that were used on the
server; they are kept as a record of how each stage was run. Colon keeps its model-selection helper scripts
under `colon/logs/`, breast cancer under `breast_cancer/pipeline/`.

### Analysis scripts (`analysis/`)

Each script starts with a comment saying what it produces and how to run it.

| Topic | Scripts |
|---|---|
| Benchmark summary and main figures | `visualisations.R`, `benchmark_summary.R`, `benchmark_figures.R`, `overall_rank.R`, `basic_metrics.R`, `compare_tissues_dissected.R`, `fetal_liver_leaders.R` |
| Are method differences meaningful (ANOVA, Hsu's MCB) | `mcb_analysis.R`, `mcb_component_table.R`, `mcb_component_figure.R`, `anova_table.R` |
| Uncertainty by subsampling | `bootstrap_scores.py`, `bootstrap_extra.py`, `bootstrap_genes.py`, `bootstrap_labels.py`, `bootstrap_bio6_summarise.py`, `bootstrap_metric_boxplots.R` |
| Additional metrics | `extra_metrics.py`, `sas_scores.py`, `new_metrics_boxplots.R` |
| Input type and orthologue removal | `orthologue_input.R`, `input_robustness.R`, `input_summary_figure.R`, `native_shared_genes.R`, `native_slopes.R` |
| Sensitivity to the labels | `label_robustness.R`, `label_sensitivity_select.R`, `label_sensitivity_scatter.R`, `label_map_alt.py`, `tau_reference.R`, `tau_extra.py` |
| Annotation quality and EBCC | `annotation_auroc_heldout.R`, `vote_vs_ebcc.R`, `ebcc_support.R`, `ebcc_support_export.py`, `seurat_vs_ebcc.R`, `seurat_confidence.R`, `best_annotator_select.R`, `best_annotator_plot.R`, `marker_auroc_figure.R`, `marker_auroc_seurat.R`, `appendix_confidence_auroc.R`, `fl_seurat_overrides.R`, `fl_seurat_overrides_flow.R`, `annotator_agreement_export.py`, `annotator_agreement_plot.R`, `singler_cross_counts.py`, `singler_cross_barplot.R` |
| Markers, panel coverage and sparsity | `marker_table.R`, `trusted_table.R`, `marker_dotplot_compute.py`, `marker_dotplot_figure.R`, `check_marker_coverage.R`, `check_genes_per_cell.R`, `check_zero_genes.py`, `check_sparsity_concentration.R`, `check_singler_genes.R`, `check_singler_handpicked.R`, `singler_genes_figure.R`, `singler_handpicked_figure.R`, `sparsity_coverage_figures.R` |
| Comparison with published benchmarks | `literature_comparison.R`, `literature_dumbbell.R` |
| Appendix figures | `dissection_figure_compute.py`, `dissection_figure.R`, `bm_umap_compute.py`, `bm_umap_figure.R`, `appendix_unintegrated_umap_compute.py`, `appendix_unintegrated_umap.R`, `appendix_pca_variance_compute.py`, `appendix_pca_variance_figure.R` |
| One-off job scripts | `run_wtau.sh`, `run_s24.sh`, `run_ngenes.sh`, `run_umap_recompute.sh` |

## What is not in this repository

- **Raw and processed data** (about 620 GB): raw spatial transcriptomics outputs, annotation references, Seurat and
  AnnData objects, integrated embeddings. The pipelines expect them under `<tissue>/data/` (raw inputs in
  `<tissue>/data/00_raw_data/`, references in `<tissue>/data/02_annotation_refs/`); the data sources are given
  in the thesis.
- **UMAP figures of every integration** (`figures/<tissue>/umap/`, about 2 GB); `analysis/visualisations.R`
  regenerates them.
- **Per-repeat intermediate files** of the bootstrap analyses (the drawn cell subsamples, and the per-repeat
  tables behind `bootstrap_extra/extra_long.csv` and `bootstrap_genes/genes_long.csv`), and all run logs.
- **Supplementary tables of the published benchmarks** read by `analysis/literature_comparison.R` from
  `scrna_benchmark/`. They belong to the publishers; download them from the papers' supplementary material
  (file names are at the top of that script).

## Running the code

The code was written for and run on a Linux server. Two settings replace the server's folder names:

```bash
export PROJECT_ROOT=/path/to/this/repository      # folder that holds library.R and the tissue folders
export CONDA_ENVS=/path/to/miniconda3/envs        # folder that holds the conda environments below
```

- Analysis scripts are run from the repository root, for example `Rscript analysis/mcb_analysis.R` or
  `python analysis/tau_extra.py`.
- Tissue scripts are run from their tissue folder, for example `cd colon && bash run_04.sh`.
- If `PROJECT_ROOT` is not set, R and Python scripts fall back to the folder they are run from (R) or to
  the repository that contains them (Python); the shell scripts stop and ask for it.
- `R_LIBS_PROJECT` can name extra R package folders to put first on the library path (optional).

### Environments

| File | Used for |
|---|---|
| `envs/R-packages.txt` | R and the packages loaded by `library.R` (stages 01-04, `visualisations.R`) |
| `envs/torch.yml` | scVI, scANVI, Scanorama, BBKNN; AnnData conversion from R |
| `envs/scib.yml` | Evaluation (stage 05) and the Python analysis scripts |
| `envs/samap.yml` | SAMap |
| `envs/tacco_env.yml` | TACCO annotation (fetal liver) |
| `envs/r45.yml` | Statistics and figure scripts that do not need Seurat (`Rscript --vanilla`) |

Create one with `conda env create -f envs/scib.yml`.

### What can be rerun from this repository alone

Scripts that read only saved score tables run without the large data, for example
`analysis/mcb_analysis.R`, `analysis/mcb_component_table.R`, `analysis/input_summary_figure.R`,
`analysis/marker_auroc_figure.R`, `analysis/best_annotator_plot.R` and `analysis/tau_extra.py`. These were
rerun from a copy of this repository and reproduced the result tables in `figures/` exactly.
Everything from stage 01 to 05, and any script that reads cell-level objects, needs the data.

## Note on paths

In the working copy on the server, the analysis scripts sat loose in the top folder and all scripts used
the server's absolute paths. For this repository they were moved into `analysis/` and the absolute paths
were replaced by `PROJECT_ROOT` / `CONDA_ENVS`. No analysis logic was changed. The full pipelines
(stages 01-05) were not rerun after this change.
