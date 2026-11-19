# Pre-processing code: verbatim excerpts

Copied straight from the scripts (path and line numbers given for each block).
Seurat functions run with their defaults wherever no argument is shown:
`NormalizeData()` = LogNormalize, scale.factor 1e4; `ScaleData()` = centre + scale, clip at 10
(with no variable features set, Seurat scales all genes; checked: every panel gene is in
`scale.data`, max 10); `RunPCA()` = 50 PCs, `seed.use = 42`; `RunUMAP()` = uwot, cosine,
30 neighbours, `min.dist = 0.3`, `seed.use = 42`.

## 1. Normalisation and scaling, per species (02 scripts)

The last call before each 02 object is saved; this is the version that goes into 03.
All panel genes are used (`VariableFeatures = rownames`); no highly variable gene selection.

### Bone marrow

`bone_marrow/02_human_annotate.R`, lines 227–230

```r
VariableFeatures(human) = rownames(human)
human = human |> 
NormalizeData()  |> 
ScaleData()
```

`bone_marrow/02_mouse_annotate.R`, lines 231–234

```r
VariableFeatures(mouse) = rownames(mouse)
mouse = mouse |> 
NormalizeData()  |> 
ScaleData()
```

### Colon

`colon/02_human_annotated.R`, lines 260–263

```r
VariableFeatures(human) = rownames(human)
human = human |>
NormalizeData()  |>
ScaleData()
```

`colon/02_mouse_annotate.qmd`, lines 516–519

```r
VariableFeatures(mouse) = rownames(mouse)
mouse = mouse |> 
NormalizeData()  |> 
ScaleData()
```

### Fetal liver

`fetal_liver/02_human_annotate.R`, lines 8–11

```r
VariableFeatures(human) = rownames(human)
human = human |>
    NormalizeData() |>
    ScaleData()
```

`fetal_liver/02_mouse_annotate.R`, lines 8–11

```r
VariableFeatures(mouse) = rownames(mouse)
mouse = mouse |>
    NormalizeData() |>
    ScaleData()
```

### Breast cancer

No `VariableFeatures()` call, so `ScaleData()` scales all genes. The PCA/UMAP/neighbours here are per species, used for annotation only (dropped in 03 by `drop_reductions()`).

`breast_cancer/02_human_annotate.qmd`, lines 51–57

```r
human = human |>
  NormalizeData() |>
  ScaleData() |>
  RunPCA(features = rownames(human)) |>
  RunUMAP(dims = 1:20, return.model = TRUE) |>
  FindNeighbors(reduction = 'pca', dims = 1:20)
SaveSeuratRds(human, 'data/02_annotated_data/02_human_annotated.rds')
```

`breast_cancer/02_mouse_annotate.qmd`, lines 58–64

```r
mouse = mouse |> 
  NormalizeData() |> 
  ScaleData() |> 
  RunPCA(features = rownames(mouse)) |> 
  RunUMAP(dims = 1:20, return.model = TRUE) |> 
  FindNeighbors(reduction = 'pca', dims = 1:20)
SaveSeuratRds(mouse, 'data/02_annotated_data/02_mouse_annotated.rds')
```

## 2. Merge, PCA and UMAP of the merged data (03 scripts)

`concatenate()` (source.R) merges the two species on one-to-one orthologues; `merge()`
keeps each species' own `data` and `scale.data`, so in bone marrow, fetal liver and breast
cancer the merged PCA runs on **per-species scaled** values. Colon rescales the merged
object jointly (`ScaleData(features = rownames(concat))`) before PCA.

### Bone marrow (full and dissected; also clusters)

`bone_marrow/03_concatenate_data.R`, lines 25–42

```r
# concatenate on one-to-one orthologues, drop empty cells, embed and (if
# resolutions are given) cluster
make_concat = function(human, mouse, resolution = NULL) {
  concat = concatenate(human, mouse, homology_tbl)
  VariableFeatures(concat) = rownames(concat)
  concat = JoinLayers(concat)
  concat$nCount_Xenium = colSums(x = concat[['Xenium']]$counts)
  concat = concat |>
    subset(nCount_Xenium > 0) |>
    RunPCA() |>
    RunUMAP(dims = 1:20, uwot.init = Sys.getenv('UMAP_INIT', 'spectral'))   # UMAP_INIT=random if spectral hangs (genes2)
  if (!is.null(resolution)) {
    concat = concat |>
      FindNeighbors(dims = 1:20) |>
      FindClusters(resolution = resolution, algorithm = 4, random.seed = 1)
  }
  concat
}
```

`bone_marrow/03_concatenate_data.R`, lines 64–73

```r
if ('full' %in% stages) {
  concat = make_concat(human, mouse, resolution = c(0.4, 0.6, 0.8, 1.2))
  save_both(concat, 'data/03_concatenated_data/', '03_concatenated_data_full')
}

# subsection
if ('subset' %in% stages) {
  concat = make_concat(human, mouse_subset, resolution = c(0.2, 0.4, 0.6, 0.8, 1.2))
  save_both(concat, 'data/03_concatenated_data/', '03_concatenated_data_subset')
}
```

### Colon (full, dissected and gene-removal levels all use `prepare_concat`)

`colon/03_concatenate.qmd`, lines 134–148

```r
# drop cells with no counts on the shared genes, then scale both species jointly
# (merge() keeps each species' own scale.data) so the baseline stays unintegrated
prepare_concat = function(concat) {
  VariableFeatures(concat) = rownames(concat)
  concat = JoinLayers(concat)
  concat$nCount_Xenium = colSums(x = concat[['Xenium']]$counts)
  concat = concat |>
    subset(nCount_Xenium > 0)
  concat |>
    ScaleData(features = rownames(concat)) |>
    RunPCA(features = rownames(concat), npcs = min(50, nrow(concat) - 1)) |>
    # PCA init: with few genes many cells are identical, the kNN graph falls apart and
    # uwot's spectral init can hang (the UMAP here is only for viewing)
    RunUMAP(dims = 1:min(20, nrow(concat) - 1), uwot.init = 'pca')
}
```

`colon/03_concatenate.qmd`, lines 158–167

```r
# human + full mouse
concat = prepare_concat(concatenate(human, mouse, homology_tbl))
save_concat(concat, 'data/03_concatenated_data/03_concatenated_data_full')
save_species_inputs(concat, human, mouse, 'data/03_concatenated_data/full_mouse/')
print(sprintf('full: %d genes x %d cells', nrow(concat), ncol(concat)))

# human + dissected mouse
concat = prepare_concat(concatenate(human, mouse_subset, homology_tbl))
save_concat(concat, 'data/03_concatenated_data/03_concatenated_data_subset')
save_species_inputs(concat, human, mouse_subset, 'data/03_concatenated_data/dissected_mouse/')
```

### Fetal liver

`fetal_liver/03_concatenate.R`, lines 52–62

```r
# full section
concat = concatenate(human, mouse, homology_tbl)
VariableFeatures(concat) = rownames(concat)
concat = JoinLayers(concat)
concat$nCount_Xenium = colSums(x = concat[['Xenium']]$counts)
concat = concat |>
  subset(nCount_Xenium > 0)

concat = concat |>
  RunPCA() |>
  RunUMAP(dims = 1:20)
```

`fetal_liver/03_concatenate.R`, lines 70–80

```r
concat = concatenate(human, mouse_subset, homology_tbl)
VariableFeatures(concat) = rownames(concat)
concat = JoinLayers(concat)
concat$nCount_Xenium = colSums(x = concat[['Xenium']]$counts)
concat = concat |>
  subset(nCount_Xenium > 0)

concat = concat |>
  RunPCA() |>
  RunUMAP(dims = 1:20)

```

### Breast cancer

`breast_cancer/03_concatenate.qmd`, lines 29–37

```r
concat = concatenate(human, mouse, homology_tbl,
                     tech_human = 'MERFISH', tech_mouse = 'MERFISH')
VariableFeatures(concat) = rownames(concat)
concat = JoinLayers(concat)
concat$nCount_MERFISH = colSums(x = concat[['MERFISH']]$counts)
concat = concat |> 
  subset(nCount_MERFISH > 0) |> 
  RunPCA(npcs = min(50, nrow(concat) - 1)) |>     # PCA + UMAP of the unintegrated data, used by the figures
  RunUMAP(dims = 1:20)
```

`breast_cancer/03_concatenate.qmd`, lines 179–187

```r
concat = concatenate(human, mouse_subset, homology_tbl,
                     tech_human = 'MERFISH', tech_mouse = 'MERFISH')
VariableFeatures(concat) = rownames(concat)
concat = JoinLayers(concat)
concat$nCount_MERFISH = colSums(x = concat[['MERFISH']]$counts)
concat = concat |> 
  subset(nCount_MERFISH > 0) |> 
  RunPCA(npcs = min(50, nrow(concat) - 1)) |>     # PCA + UMAP of the unintegrated data, used by the figures
  RunUMAP(dims = 1:20)
```

## 3. Gene-removal datasets (genes2/4/6/8): renormalise per species, then merge as above

Counts unchanged; each species is log-normalised and scaled again on the genes that remain.


`bone_marrow/03_concatenate_data.R`, lines 89–114

```r
# counts are unchanged; data/scale.data must be recomputed on the genes that remain
renormalize = function(obj) {
  obj$nCount_Xenium = colSums(x = obj[['Xenium']]$counts)
  obj = obj |>
    subset(nCount_Xenium > 0)
  VariableFeatures(obj) = rownames(obj)
  obj |>
    NormalizeData() |>
    ScaleData()
}

for (p in as.numeric(gene_stages)) {
  out_dir = paste0('data/genes', p * 10, '/')
  removed = shared_genes[removal_order[seq_len(round(p * nrow(shared_genes)))], ]
  human_p = human[!rownames(human) %in% removed$hsapiens_homolog_ensembl_gene, ] |>
    renormalize()
  mouse_p = mouse_subset[!rownames(mouse_subset) %in% removed$mmusculus_homolog_ensembl_gene, ] |>
    renormalize()
  concat_p = make_concat(human_p, mouse_p)
  print(str_glue('{out_dir}: removed {nrow(removed)} of {nrow(shared_genes)} shared genes, ',
                 '{nrow(concat_p)} remain'))

  save_both(human_p, out_dir, '02_human_annotated')
  save_both(mouse_p, out_dir, '02_mouse_annotated')
  save_both(concat_p, out_dir, '03_concatenated_data')
}
```

`colon/03_concatenate.qmd`, lines 184–193

```r
# levels run one after another: forking them (mclapply) segfaults in RunPCA / RunUMAP
for (frac in c(0.2, 0.4, 0.6, 0.8)) {
  drop = gene_order[seq_len(round(frac * n_shared))]
  # each species renormalised (and rescaled, as in 02) on its remaining genes, then concatenated as above
  renorm = function(x) x |> NormalizeData() |> ScaleData(features = rownames(x))
  human_p = renorm(human[!rownames(human) %in% drop, ])
  mouse_p = renorm(mouse_subset[!rownames(mouse_subset) %in% hs_to_mm[drop], ])
  stopifnot(nrow(human) - nrow(human_p) == length(drop),
            nrow(mouse_subset) - nrow(mouse_p) == length(drop))
  concat_p = prepare_concat(concatenate(human_p, mouse_p, homology_tbl))
```

`fetal_liver/03_concatenate.R`, lines 102–127

```r
done = lapply(p, FUN = function(frac) {
    drop = gene_order[seq_len(round(frac * n_shared))]
    # renormalise each species on its remaining genes (log-normalisation and scaling
    # depend on which genes are present), then concatenate them as for full / dissected
    renormalise = function(obj) {
      VariableFeatures(obj) = rownames(obj)
      obj |> NormalizeData(verbose = FALSE) |> ScaleData(verbose = FALSE)
    }
    human_p = renormalise(human[!rownames(human) %in% drop, ])
    mouse_p = renormalise(mouse[!rownames(mouse) %in% hs_to_mm[drop], ])
    stopifnot(nrow(human) - nrow(human_p) == length(drop),
              nrow(mouse) - nrow(mouse_p) == length(drop))
    concat_p = concatenate(human_p, mouse_p, homology_tbl)
    stopifnot(setequal(rownames(concat_p), setdiff(rownames(concat), drop)))
    VariableFeatures(concat_p) = rownames(concat_p)
    concat_p = JoinLayers(concat_p)
    # cells left with no counts on the kept shared genes are removed everywhere
    concat_p$nCount_Xenium = colSums(x = concat_p[['Xenium']]$counts)
    concat_p = concat_p |>
      subset(nCount_Xenium > 0)
    concat_p = concat_p |>
      RunPCA(npcs = min(50, nrow(concat_p) - 1)) |>
      # PCA init: with few genes many cells are identical, the kNN graph falls apart and
      # uwot's spectral init hangs or segfaults in RSpectra (UMAP here is only for viewing)
      RunUMAP(dims = 1:min(20, nrow(concat_p) - 1), uwot.init = 'pca')

```

`breast_cancer/03_concatenate.qmd`, lines 260–285

```r
# counts are unchanged; data/scale.data are recomputed on the genes that remain
renormalize = function(obj) {
  obj$nCount_MERFISH = colSums(x = obj[['MERFISH']]$counts)
  obj = obj |>
    subset(nCount_MERFISH > 0)
  VariableFeatures(obj) = rownames(obj)
  obj |>
    NormalizeData() |>
    ScaleData()
}

# species objects carry the hsapiens_/mmusculus_ prefixes; concatenate() adds them itself
make_concat = function(human, mouse) {
  human = RenameCells(human, new.names = sub('^hsapiens_', '', colnames(human)))
  mouse = RenameCells(mouse, new.names = sub('^mmusculus_', '', colnames(mouse)))
  concat = concatenate(human, mouse, homology_tbl, tech_human = 'MERFISH', tech_mouse = 'MERFISH')
  VariableFeatures(concat) = rownames(concat)
  concat = JoinLayers(concat)
  concat$nCount_MERFISH = colSums(x = concat[['MERFISH']]$counts)
  concat = concat |>
    subset(nCount_MERFISH > 0)
  # PCA + UMAP of the unintegrated data, used by the figures
  concat |>
    RunPCA(npcs = min(50, nrow(concat) - 1)) |>
    RunUMAP(dims = 1:20)
}
```

## 4. Related: joint PCA used as the unintegrated baseline in evaluation (05)

Not part of 03, but it is the PCA that the PCR_batch baseline uses, scaled jointly over
both species (all tissues have the same function in their utils.py).

`bone_marrow/utils.py`, lines 22–31

```python
# PCA of the unintegrated data: log-normalised `data` layer, scaled jointly over both
# species. Not the `scale.data` layer: it was scaled within each species in 02, which
# removes most of the species effect (PCR_batch baseline ~0, so most methods scored 0).
# X is raw counts, and the stored X_pca was computed before genes were dropped.
def joint_pca(adata):
    data_layer = next((l for l in adata.layers if l.startswith('data')), None)
    if data_layer is not None:
        adata.X = adata.layers[data_layer].copy()
    sc.pp.scale(adata, max_value = 10)
    sc.pp.pca(adata)
```
