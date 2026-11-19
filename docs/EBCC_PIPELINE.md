# EBCC model selection, downstream rerun and figures

Reference procedure, established on `bone_marrow`. Apply it the same way in the
other tissues (colon, breast_cancer, breast_cancer_cross, fetal_liver), mapping
the steps onto each project's own files. bone_marrow paths are given as examples.

## 0. Ground rules
- **M is chosen by held-out log-likelihood** (section 3, joint (M, frac) grid).
  **frac is then chosen by marker AUROC with M fixed** (section 3b), using a 50/50
  cell split: pick on one half, report on the held-out half.
- Why two criteria: the held-out log-likelihood measures agreement with the
  annotators, so it always favours a weak trusted prior (it picked frac = 0 for
  both bone_marrow species, which switches the `trusted` lists off and lowered
  marker agreement). frac is the knob that encodes marker-based trust, so it is
  chosen by marker agreement instead.
- State this plainly when reporting: frac was tuned for marker agreement, and the
  same markers informed the `trusted` lists, so the marker AUROC is not an
  independent validation of EBCC.
- Do not change the `trusted` lists or the marker lists during this procedure;
  only M and frac are selected.
- Never edit a shell script while bash is running it (bash reads scripts
  incrementally; inserted bytes corrupt the run). Patch after it exits.
- Run long jobs in `screen`. The server is shared: check `uptime` and cap threads.

## 1. Prerequisites in each `02_<species>_annotate.R` (check before starting)
- Reference `X` from h5ad is usually already log-normalised: use
  `counts = expm1(X)`, drop `ADT-` features, then `NormalizeData()`; no `ScaleData()`
  on references.
- Both references must share one label set (bone_marrow: Mast -> Granulocyte,
  Mac -> Monocyte).
- Cross-species reference: `concatenate(ref, query, homology_tbl)` with
  `tech_mouse` matching the mouse object's assay (Xenium query, not RNA).
- `ann_cols = c("celltype_seurat", "celltype_singleR", "celltype_seurat_cross",
  "celltype_singleR_cross")` (first two: same-species reference; `_cross`: other species').
- `classes` built from **all four** annotator columns.

## 2. Held-out log-likelihood (in the 02 scripts)
- `set.seed(1)`; sample 10% of cells; for each, hide one randomly chosen
  annotator's label (`meta_ho`, `ho`). The same held-out set is used for every fit.
- Prior: `make_beta(M)` = `make_beta_prior(classes, ann_cols)` (diag 4, off 1),
  then for each trusted (annotator, class) **diag = 4 + frac * n_k / M**, where
  n_k = number of cells that annotator gave that label. Skip trusted classes the
  annotator never called.
- `heldout_ll(M, R = 3)`: `run_ebcc(meta_ho, ..., keep_rho = TRUE)`; per hidden
  label, `log( sum over class-subtypes of P(cell in subtype) x P(annotator gives that label) )`.

## 3. Joint (M, frac) grid
- Grid: M in {1, 2, 3, 5, 10} x frac in {0, 0.25, 0.5, 1, 2}; 3 restarts per fit;
  per species. Run one process per (species, frac) in parallel, each over all M
  (OMP/OPENBLAS/N_THREADS = 4 with 10 jobs).
- Rule: per (M, frac) mean and SE (sd/sqrt(n)) of the held-out log-likelihoods.
  Best = highest mean. Among combinations within 1 SE of the best, choose the
  **smallest M, then the smallest frac**.
- Save:
  - `data/02_annotated_data/<species>_ebcc_heldout_grid.rds`: all per-cell held-out values, summary, rule text
  - `data/02_annotated_data/<species>_ebcc_heldout_grid.csv`: M, frac, mean_ll, se, n_heldout, within_1se_of_best, best, chosen
  - `data/02_annotated_data/<species>_ebcc_heldout_grid_chosen.csv`: species, M, frac, changed, previous_M, previous_frac
  - `figures/<tissue>/annotation/ebcc_heldout_grid_<species>.png`: heatmap (thin outline = within 1 SE, thick = chosen)
- **Take only M from this grid** (the chosen row's M); frac is chosen in 3b.

## 3b. frac by marker AUROC (M fixed)
- With M fixed from section 3, run a **full** EBCC fit (`R = 10`, all cells) for
  each frac in {0, 0.25, 0.5, 1, 2, 5}, per species, in parallel. Save each fit
  (label, posterior, par, elbo_restarts) so the chosen one becomes the final fit
  without refitting.
- Marker score per cell and type = mean log-normalised expression of that type's
  panel markers (markers from the 02 script's `markers = list(...)`).
  AUROC per type = P(cell labelled with the type scores higher than a cell
  labelled otherwise). Summary = mean over the types the labels call.
- `set.seed(1)`; split cells 50/50 at random. **Choose frac with the highest mean
  AUROC on the selection half** (ties -> smaller frac); report the AUROC on the
  held-out half next to it.
- Save:
  - `data/02_annotated_data/<species>_ebcc_frac_auroc.csv`: frac, mean AUROC per half, number of types, chosen
  - `data/02_annotated_data/<species>_ebcc_frac_auroc_per_class.csv`: frac, half, cell type, n cells, AUROC
  - `data/02_annotated_data/<species>_ebcc_frac_auroc_chosen.csv`: species, M, frac, changed, previous_frac, rule
  - `figures/<tissue>/annotation/ebcc_frac_auroc_<species>.png`: mean AUROC vs frac, both halves, chosen marked
- Expect AUROC to rise with frac (EBCC follows the trusted annotators more); if
  the largest frac wins, say so — the choice then sits at the edge of the grid.
- "Changed" = chosen M or frac differs from that species' 02 script. If neither
  species changed, stop here.

## 4. Final fit (each changed species)
- Use the saved fit for the chosen (M, frac) from 3b (no refit needed); it was
  `run_ebcc(<obj>@meta.data, ann_cols, M = M, R = 10, beta = make_beta(M))`.
- Save `ebcc_label`, `ebcc_maxpost`, `ebcc_entropy`, `celltype = ebcc_label`;
  `<species>_ebcc_params.rds`, `<species>_ebcc_posterior.rds`,
  `02_<species>_annotated.rds/.h5ad`; EBCC marker dotplot to
  `figures/<tissue>/annotation/<species>_dotplot_ebcc.png`.
- Write the chosen `frac = ...` and `M = ...` (or the hard-coded M in the
  `run_ebcc` call) back into the 02 script, with a comment pointing to the
  selection CSV.

### bone_marrow outcome (for reference)
- Held-out log-likelihood grid: human best M=10/frac=0, chosen M=5, frac=0;
  mouse best M=10/frac=0, chosen M=3, frac=0. Log-likelihood fell steadily with
  frac at every M.
- With frac = 0, EBCC marker AUROC dropped (human 0.72 -> 0.70, mouse 0.75 -> 0.73),
  which is why frac is now chosen by marker AUROC (3b) with M = 5 (human), 3 (mouse).

## 5. Downstream (only if something changed)
1. Mouse changed: redo the mouse dissection (the tissue-matched cut from 5b) ->
   `02_mouse_subset_annotated.rds/.h5ad`.
2. Fingerprint each concatenated dataset (md5 of cell and gene names), then rerun
   all of 03 (`./run_03_04.sh 03`).
3. Cells/genes unchanged: rerun **only scANVI (04.08)**, the only label-dependent
   method (it loads the saved scVI model `04.07_scVI_model`). Otherwise rerun all
   of 04 (`./run_03_04.sh 04`).
4. Rerun all of 05 (`./run_05.sh`).
5. Delete this tissue's UMAP and spatial figure folders and rerun
   `Rscript visualisations.R <tissue>` (tissue names on the command line run only
   those; never redraw other tissues' figures while their sessions are writing data).

## 5b. Dissection: match the mouse region to the human tissue
The `dissected` dataset (and the genes2/4/6/8 levels built from it) should hold
only the mouse compartment that matches what the human sample actually is. A cut
that still contains tissue the human sample lacks (cartilage, bone, capsule,
muscle, other organs) gives the integration mouse-only cell populations and makes
the cross-species comparison unfair.
1. **Decide what the human sample is** (e.g. adult marrow aspirate/core, colon
   mucosa, tumour core, fetal liver parenchyma) and which mouse compartment
   matches it. Write this down in one sentence; it is the justification.
2. **Map the mouse section**: bin `GetTissueCoordinates()` into ~250 x 250
   squares. Per bin, take the fraction of cells with target-compartment EBCC labels
   and the mean expression of off-target marker sets, which do not depend on the
   labels. bone_marrow used erythroid/granulocyte/monocyte/MkP labels as
   blood-forming, cartilage `Col2a1/Acan/Comp/Ucma` and bone `Col1a1/Col1a2/Runx2`;
   marrow bin = >= 70% blood-forming and mean cartilage < 1. Plot the bin map.
3. **Choose the region from the map**: a simple, reproducible shape (an
   axis-aligned rectangle in `GetTissueCoordinates` units; note that `y` is the
   horizontal axis and `x` the vertical axis in `ImageDimPlot`). It should keep
   the target compartment and leave out the off-target tissue. Prefer one
   contiguous anatomical unit over a bin mask (easier to defend and to redo).
4. **Report before / after composition** (EBCC label %, whole section vs cut, plus
   the human composition) to `data/02_annotated_data/mouse_dissection_composition.csv`,
   and save the cut over the whole section to
   `figures/<tissue>/spatial/mouse_dissection.png`. Off-target labels (e.g. MSC,
   or immune labels that sit on cartilage) should drop. Any composition difference
   that remains and is real biology (bone_marrow: neonatal mouse marrow is
   erythroid-rich) is stated as a limitation and is **not** trimmed away.
5. Put the cut in the notebook/script that writes `02_mouse_subset_annotated.*`
   (the coordinates are written in the code, not drawn by hand), and export the
   h5ad with `RETICULATE_PYTHON` set to the `torch` env (`as.anndata` needs
   anndata). Export the same variable before running 03 from a screen or chain
   script, or 03 fails at its h5ad export.
6. Rerun downstream for **`dissected genes2 genes4 genes6 genes8` only**: the full
   dataset does not change. Both runners take a `DATASETS` override:
   `Rscript 03_concatenate_data.R inputs` (rewrites `03_mouse_subset.h5ad`, the
   SAMap input; skipping it leaves SAMap on the old cells with no error), then
   `subset`, then the gene levels, then
   `DATASETS="dissected genes2 genes4 genes6 genes8" ./run_03_04.sh 04`, the same
   for `./run_05.sh`, then the 9b sensitivity run (if this tissue has one) with the same
   `DATASETS`, then delete `figures/<tissue>/umap/{dissected,genes2,genes4,genes6,genes8}`
   (`visualisations.R` skips figures that already exist) and run
   `Rscript visualisations.R <tissue>`. All of 04 reruns because
   the cells change.
   (Check first that your tissue's gene levels are built from the dissected
   subset, as in bone_marrow; if they come from the full data, rerun only
   `dissected`.)
- bone_marrow (P8 mouse, both femurs vs adult human marrow): the old cut
  (`1700 <= y <= 4700`) kept epiphyses/cartilage and cortical stroma (MSC 21%,
  T/NK on cartilage). New cut = right femur marrow cavity,
  `2500 <= y <= 3500 & 500 <= x <= 4750`, 24,815 / 52,045 cells. MSC 21 -> 7%,
  T cell 2.2 -> 0.1%, NK/ILC 3.3 -> 0.5%, erythroid 46 -> 65% (neonatal biology,
  reported as a limitation).

## 6. Pipeline conventions
- Datasets: full, dissected, genes2/4/6/8 (20/40/60/80% of shared one-to-one
  orthologues removed; nested removal; per-species data renormalised on the
  remaining genes). Listed once in `datasets.py` / `integration_sets.R`; every
  04/05 script takes dataset names as arguments (default: all).
- `run_03_04.sh` must export `N_THREADS` equal to the OpenMP/BLAS cap
  (`default_threads()` in `source.R` reads it; Harmony refuses a mismatch).
- 04.11 (LIGER UINMF) must finish before 04.12-04.17 and 05 (`uinmf_H_raw.*`,
  `eval_cells.txt`).
- `eval_cells.txt` = cells in the concatenated data whose UINMF H row has
  **>= 2 non-zero factors** (RPCA UINMF gives NA otherwise). Every method and both
  baselines are scored on exactly these cells (`utils.align(..., eval_cells(...))`).
- Python envs: `torch` (scVI, scANVI, Scanorama, BBKNN), `samap` (SAMap, CUDA
  backend via CuPy), `scib` (05).

## 6b. Methods on an scVI latent (04.18-04.24, 05.18-05.24; added 2026-10-26)
Same idea as the UINMF variants (04.12-04.17, which run on the UINMF H matrix): run
FastMNN, CCA, Harmony, Scanorama, RPCA and BBKNN on an **scVI latent trained without
the species key** instead of PCA. It is a denoised, non-linear counterpart of PCA in
which the species are still separate, so the question is whether a VAE latent is a
better basis than PCA (or UINMF H). Do NOT use the integrated `X_scvi` from 04.07: the
methods would then correct an already-corrected space.
- Copy from `bone_marrow/`: `04.18_scVI_latent.py`, `04.19_integrate_fastMNN_scvi.R`,
  `04.20_integrate_CCA_scvi.R`, `04.21_integrate_Harmony_scvi.R`,
  `04.22_integrate_Scanorama_scvi.py`, `04.23_integrate_RPCA_scvi.R`,
  `04.24_integrate_BBKNN_scvi.py`, `05.18`-`05.24_evaluate_*_scvi*.py`, and the
  `04scvi` part of `run_03_04.sh` (plus the GPU-lane / phase-3 lines, and
  `xargs -r` in `pool()`). `source.R` (shared) already has `read_scvi_latent()` and
  `rename_reduction()`; `visualisations.R` (shared) already lists the 7 new methods
  and skips any whose 04 file a tissue does not have.
- 04.18: `scvi.model.SCVI.setup_anndata(adata, layer = 'counts')` with **no
  batch_key**, `n_latent = 20` (matches the 20 UINMF factors), `scvi.settings.seed = 1`;
  writes `scvi_latent.csv` (cells x 20, concatenated cell names with species prefix)
  and `04.18_scVI_latent.h5ad` (`X_scvi_latent`, `X_umap_scvi_latent`).
- 04.19-04.24 pass the latent to the existing `integrate_*` helpers through their
  `H` argument (`read_scvi_latent()` returns the same per-species list as the UINMF H),
  with `W = matrix(numeric(0), 0, 0)`, then rename the UINMF-named reductions
  (`mnn_scvi`, `cca_scvi`, `rpca_scvi`; Harmony stays `harmony`).
- 05.18 scores the latent itself as "Unintegrated scVI-latent" (the baseline for the
  family); all are scored on `eval_cells.txt` like every other method.
- Method names end in ` scVI-latent`; in the score figures they share the base
  algorithm's colour and get their own shape (open square).
- Run: `DATASETS="..." ./run_03_04.sh 04scvi` (04.18 on the GPU, then 04.19-04.24),
  then `SCRIPTS="05.1[89]*.py 05.2*.py" ./run_05.sh` (only the 7 new scores; copy the
  `SCRIPTS` line from `bone_marrow/run_05.sh`), and the 9b sensitivity run with the same
  `SCRIPTS` if the tissue has one, then `Rscript visualisations.R <tissue>`. scVI training is ~8 s/epoch
  on a busy machine: put datasets on both GPUs if free.
- Report these as a separate "scVI basis" comparison alongside the main ranking
  (the main table had 19 methods; with these it has 26).

## 7. Evaluation (05, `utils.py`)
- scib metrics, `kBET_ = False` (no separate kBET script).
- Batch = mean of `PCR_batch`, `ASW_label/batch`, `iLISI`, `graph_conn`;
  bio = mean of NMI, ARI, `ASW_label`, isolated-label F1/silhouette, `cLISI`;
  missing values ignored.
- `LABEL = 'ebcc_label'`, `BATCH = 'orig.ident'`, always from the 03 data.
- **Graph-only methods scored on the wrong graph (bug fixed 2026-10-26).** For
  `output_type = 'knn'` (SAMap, BBKNN, BBKNN UINMF, BBKNN scVI-latent) `score()` passed
  `embed = 'X_pca'`. scib's clustering then compares that with
  `uns['neighbors']['params']['use_rep']` and, if different, **rebuilds the kNN graph
  in place from `X_pca`**; NMI, ARI, graph_conn, iLISI and cLISI were then all computed on
  that rebuilt graph, not the method's. Affected: SAMap (graph built on `X`; 05.02 adds
  a PCA of SAMap's X) and BBKNN UINMF (`X_uinmf_H`; the carried-over `X_pca` is the
  unintegrated PCA, hence iLISI = 0). Plain BBKNN (built on `X_pca`) was fine.
  Check: `grep "Recompute neighbors" logs/05/*.log`.
  Fix (`utils.py`, in `score()` before `scib.me.metrics`): for knn, set
  `embed = adata_int.uns['neighbors']['params']['use_rep']`, and sort the CSR indices of
  `obsp['connectivities']` / `['distances']` (`align()` reorders cells; the graph-LISI
  C++ helper then prints "malformed matrix format" and crashes once the graph is
  actually used). Then rerun 05.02 and 05.17 (and 05.24) for every dataset, the 9b
  sensitivity run for those scripts, and the figures. bone_marrow genes8 SAMap: iLISI
  0 -> 0.23, ARI 0.28 -> 0.36; BBKNN UINMF: iLISI 0 -> 0.68.
- **PCR_batch baseline (bug fixed 2026-10-24).** The 03 `scale.data` layer is
  scaled *within each species* (inherited from 02 through `concatenate()`), which
  removes most of the species effect: in bone_marrow, species explained 0.15-0.7%
  of baseline variance vs 5-7% with joint scaling. PCR_batch =
  (before - after) / before, clipped at 0, so with that baseline most methods
  scored exactly 0 and a few got high scores, which swung the batch ranking.
  Fix (`utils.py`): `joint_pca(adata)` = `data` layer -> `sc.pp.scale(max_value=10)`
  across both species -> `sc.pp.pca`; `score()` calls it on `raw`, and
  `05.00_evaluate_raw.py` uses it for the Unintegrated embedding (instead of
  `sc.pp.pca(layer='scale.data_Xenium')`). Then rerun all 05 (and the 9b
  sensitivity run) and the score figures.
  Check: count PCR_batch == 0 over the 05 CSVs before and after; with the fix
  only methods whose embedding is genuinely more species-associated than the raw
  data stay at 0 (e.g. scANVI on bone_marrow dissected; PCR is known to penalise
  embeddings organised by cell type when cell-type composition differs between
  species).
- SAMap outputs have no species prefix: add `hsapiens_`/`mmusculus_` before aligning.
- Graph-only methods (BBKNN, BBKNN UINMF, SAMap) have no ASW/PCR: their means use
  fewer metrics. scANVI is trained on `ebcc_label` (label-informed).

## 8. Figures (`visualisations.R`, one config entry per tissue)
All under `figures/<tissue>/{scores,umap,spatial,annotation}/`, never loose in
the project folder (including the 02 scripts' per-annotator dotplots).
- Score scatter: x = batch, y = bio; colour = base algorithm (method minus
  " UINMF"); shape = UINMF or not; dashed lines at **plain Unintegrated**. One PNG
  per case + one combined; scores also saved as CSV.
- Pareto: on the front if no other method in that case is >= on both scores and
  > on one; staircase `geom_step(direction = 'hv')`; front points labelled, rest dimmed.
- UMAP: each method's **own** obsm key read from the 04 h5ad (e.g. `X_umap_mnn`,
  not the carried-over `X_umap`); only `eval_cells.txt` cells; zoom to the central
  99% (count outside in the subtitle); species | cell type side by side.
  Unintegrated UINMF: uwot on the H matrix (cosine, 30 neighbours).
- Spatial: `ImageDimPlot` per annotator; do **not** pass `border.color = NA`
  (hides centroids); one colour per label across annotators; NA -> grey "Unassigned".
- Marker AUROC: per cell type, AUROC of mean marker expression for cells with
  that label vs all others; markers read from the 02 script's `markers = list(...)`;
  one PNG per species + CSV; summary row "Mean (all called types)".

## 8b. Runtime (approximate, from the run logs; added 2026-10-26)
Runtime is reported as a rough guide in tiers (minutes vs tens of minutes), taken from
existing logs, with no new runs. It cannot rank methods within a tier.
- Copy `bone_marrow/08_runtime_from_logs.py` and point it at this tissue's logs: the
  last complete run of each dataset (a `run_03_04.sh`-style log with
  `HH:MM:SS start 04.xx_<dataset>` / `done` lines). Use the latest run on the current
  data for each dataset; reruns with refreshed inputs (e.g. SAMap after the dissection
  fix) supersede older ones. Jobs logged without timestamps: use the job log's file
  creation -> last-write time (`stat -c %W` / mtime), and skip any log that was
  rewritten by a second run (its creation time is from the first run).
- Output: `data/08_timing/runtime_from_logs.csv`, one row per method x dataset with
  columns `script, method, input (PCA / UINMF / scVI-latent / own), device (cpu/gpu),
  needs (the one-off step a variant depends on: 04.11 for UINMF variants, 04.18 for
  scVI-latent, 04.07 for scANVI), dataset, n_cells, minutes, umap_init_minutes,
  minutes_excl_umap_init, source, note`. Same column names in every tissue so the
  tissues can be plotted together.
- **Use `minutes_excl_umap_init`.** uwot's spectral initialisation ("Initializing from
  normalized Laplacian") can stall for hours on some datasets (bone_marrow: 43 min to
  ~4 h on genes2/genes4). The UMAP is for plotting only, so that time is removed by
  reading the timestamps around that line in each R job log.
- Caveats to state in the caption: wall-clock on a shared server, including reading
  inputs, UMAP and saving outputs; up to 6 CPU + 2 GPU jobs in parallel. The full cost of
  a UINMF / scVI-latent variant is its own row + the `needs` row.
- bone_marrow tiers (dissected, ~50k cells): PCA- and UINMF-based methods 1-3 min,
  scVI-latent variants 4-5 min (+27 min latent training), SAMap 8 min (GPU), LIGER UINMF
  7 min, scVI 19 min (GPU). Most roughly double from dissected to full; LIGER UINMF
  instead grows as shared genes are removed (7 -> 21 min).
- A controlled benchmark (`08_timing.R`, `08_timing.py`, `run_timing.sh`: integration
  step only, sequential, fixed threads, 3 shuffled repeats) exists in bone_marrow but was
  not run; use it only if runtime becomes a main argument.

## 9b. Label choice and sensitivity analysis
- **EBCC labels are the primary labels** for every tissue: pre-specified, the same
  method for both species and all tissues, one harmonised label set across species.
  Do not switch a tissue/species to a single annotator because it scores higher on
  the marker AUROC: that choice would be made with the same markers that set the
  trusted lists and frac (circular), and it would make tissues inconsistent.
- Report the marker AUROC heatmap descriptively (e.g. "EBCC comparable to or above
  single annotators in mouse, below same-species Seurat in human"), not as the
  reason for the label choice.
- **Targeted robustness check** (run only where EBCC is not the best annotator):
  1. From `figures/<tissue>/annotation/marker_auroc.csv`, take each annotator's
     mean AUROC over the cell types it calls ("Mean (all called types)").
  2. Per species: if EBCC has the highest mean, keep EBCC for that species;
     otherwise use the annotator with the highest mean (its 03 obs column, e.g.
     `celltype_seurat`). If EBCC is best in both species, skip the check for
     that tissue.
  3. Rescore 05 once with that mixed label set. Only 05 is rerun (integrations do
     not use labels; scANVI stays trained on EBCC and is noted as label-informed).
     - `utils.py` reads `LABEL_SET` and `LABEL_COLS`
       (`LABEL_COLS='hsapiens=<col>,mmusculus=<col>'`); unset = EBCC, unchanged.
     - `datasets.py` sends non-default runs to a `label_<set>/` subfolder of each
       scores folder; `run_05.sh` logs to `logs/05_<set>/`.
     - Run: `LABEL_SET=best_annotator LABEL_COLS=... JOBS=6 THREADS=8 ./run_05.sh`.
       Patch these files only when no 05 run is in progress.
  4. Per dataset, Kendall tau between method rankings under EBCC and under the
     best-annotator labels (overall = 0.6 bio + 0.4 batch, and bio), plus the
     top-3 methods under each. Save
     `data/05_evaluated_score/label_sensitivity.csv`,
     `data/05_evaluated_score/label_sensitivity_scores.csv`,
     `data/05_evaluated_score/label_sensitivity_choice.json` (AUROC means and the
     choice) and `figures/<tissue>/scores/label_sensitivity.png`.
  - Framing: the strongest single annotator is used as a deliberately tough
    comparison; using the AUROC to pick it is fine because this is a robustness
    check, not how the primary labels were chosen.
  - bone_marrow: human -> same-species Seurat (0.81 vs EBCC 0.74), mouse -> EBCC
    (0.76, best).

## 9. Tissue-specific checks (do not assume they carry over)
In bone_marrow the mouse T/B/NK-ILC EBCC labels were mostly chondrocytes (no
reference class covers them), and human EBCC scored below same-species Seurat on
marker AUROC. Check each tissue with the EBCC dotplots, per-label spatial plots
and the marker AUROC before trusting its labels.
