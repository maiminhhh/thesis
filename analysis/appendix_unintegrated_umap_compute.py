# UMAPs of the unintegrated data (joint log-normalised PCA, dissected mouse case) for analysis/appendix_unintegrated_umap.R,
# recomputed with exactly the settings of the bone-marrow method figure (imported from analysis/bm_umap_compute.py: kNN 15,
# euclidean, all PCs, min_dist 0.5, spread 1, seed 0, random initial layout), instead of the Seurat RunUMAP saved in the h5ad.
# Nothing is integrated or refitted; X_pca is read as saved. UMAP is run on all cells of each file; the figure keeps the
# evaluation cells. Bone marrow uses the same input as bm_umap/Unintegrated.csv, so the two must agree (checked below).
# Usage: python analysis/appendix_unintegrated_umap_compute.py
#   -> figures/benchmark_summary/unintegrated_umap_data/<tissue>.csv (cell, species, celltype, umap_1, umap_2), settings.txt
import os, h5py, numpy as np, pandas as pd
import bm_umap_compute as bu

# ---- config ----
cases = {  # tissue: concatenated h5ad (as analysis/appendix_unintegrated_umap.R)
    'colon':         'colon/data/03_concatenated_data/03_concatenated_data_subset.h5ad',
    'bone_marrow':   'bone_marrow/data/03_concatenated_data/03_concatenated_data_subset.h5ad',
    'fetal_liver':   'fetal_liver/data/03_concatenated_data/03_concatenated_data_subset.h5ad',
    'breast_cancer': 'breast_cancer/data/dissected_mouse/03_concatenated_data.h5ad',
}
rep = 'X_pca'
out = 'figures/benchmark_summary/unintegrated_umap_data'

def col(o, k):
    d = o[k]
    return d['categories'][:].astype(str)[d['codes'][:]] if isinstance(d, h5py.Group) else d[:].astype(str)

os.makedirs(out, exist_ok=True)
for tissue, f in cases.items():
    with h5py.File(f, 'r') as h:
        o = h['obs']
        cells, species, celltype = bu.cells_of(h), col(o, 'orig.ident'), col(o, 'celltype')
        emb = h['obsm/' + rep][:]
    u = bu.umap_of(cells, emb=emb, rep=rep)
    pd.DataFrame({'cell': cells, 'species': species, 'celltype': celltype, 'umap_1': u[:, 0], 'umap_2': u[:, 1]}).to_csv(
        os.path.join(out, tissue + '.csv'), index=False)
    print(f'{tissue}: {len(cells)} cells, {rep} {emb.shape[1]} dims', flush=True)

bm_new = pd.read_csv(os.path.join(out, 'bone_marrow.csv'))
bm_fig = pd.read_csv(os.path.join(bu.out, 'Unintegrated.csv'))
print('bone marrow identical to bm_umap/Unintegrated.csv:',
      bool((bm_new.cell == bm_fig.cell).all() and np.allclose(bm_new[['umap_1', 'umap_2']], bm_fig[['umap_1', 'umap_2']])))
with open(os.path.join(out, 'settings.txt'), 'w') as fh:
    fh.write(bu.settings_text().replace('; graph-based methods: their saved connectivities', '') + f' Input: {rep} (all PCs).\n')
print('done')
