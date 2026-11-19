# Inputs for analysis/ebcc_support.R read from the 02 h5ad instead of the Seurat .rds (Seurat no longer loads under the r45 env).
# Writes, for one dataset:
#   <stem>_obs.csv.gz    cell, the four annotator columns, ebcc_label, ebcc_maxpost   (h5ad cell order = Seurat order)
#   <stem>_score.csv.gz  one column per marker group with >= 1 panel gene: mean log-normalised expression of its genes
#                        per cell (as Matrix::colMeans(X[genes, ]) on the Seurat 'data' layer in analysis/ebcc_support.R)
# Usage: python -I analysis/ebcc_support_export.py <h5ad> <markers.json> <out stem>
import sys, json
import numpy as np
import pandas as pd
import anndata as ad

h5, marker_file, stem = sys.argv[1:4]
a = ad.read_h5ad(h5)
layer = next(l for l in a.layers.keys() if l == 'data' or l.startswith('data_'))   # log-normalised, not scale.data
genes = pd.Index(a.var['gene_name'] if 'gene_name' in a.var.columns else a.var_names).astype(str)
first = {}
for i, g in enumerate(genes):
    first.setdefault(g, i)                                   # duplicated names: R's X[g, ] takes the first
X = a.layers[layer].tocsc()

markers = json.load(open(marker_file))
score = {}
for k, gs in markers.items():                                # same order as annotation_markers(a)
    cols = [first[g] for g in dict.fromkeys(gs) if g in first]   # intersect(markers, rownames(X)): unique, marker order
    if cols:
        score[k] = np.asarray(X[:, cols].mean(axis=1)).ravel()
pd.DataFrame(score).to_csv(stem + '_score.csv.gz', index=False, float_format='%.17g')

cols = ['celltype_seurat', 'celltype_singleR', 'celltype_seurat_cross', 'celltype_singleR_cross', 'ebcc_label', 'ebcc_maxpost']
obs = a.obs[cols].copy()
for c in cols[:-1]:
    obs[c] = obs[c].astype(object).where(obs[c].notna(), None)
obs.insert(0, 'cell', a.obs_names.astype(str))
obs.to_csv(stem + '_obs.csv.gz', index=False)
print(f'{stem}: {a.n_obs} cells, layer {layer}, {len(score)}/{len(markers)} marker groups scored')
