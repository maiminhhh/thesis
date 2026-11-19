# Appendix: explained variance of each species' own PCA, per tissue (dissected mouse case, evaluation cells).
# Steps as extra_metrics.species_pca / utils.joint_pca, but on one species at a time: log-normalised layer
# ('data*'), sc.pp.scale(max_value=10) within the species, PCA. All components are kept (full SVD), so the
# cumulative curve reaches 100%; the first 50 ratios equal those of sc.pp.pca's default 50 components.
# Nothing is integrated or rescored.
# Usage: python analysis/appendix_pca_variance_compute.py
#   -> figures/benchmark_summary/pca_variance.csv (tissue, species, pc, var_ratio, cum_ratio, n_cells, n_genes)
import numpy as np, pandas as pd, scanpy as sc
from sklearn.decomposition import PCA

# ---- config ----
cases = {  # tissue: (concatenated h5ad, directory with eval_cells.txt), as analysis/appendix_unintegrated_umap.R
    'colon':         ('colon/data/03_concatenated_data/03_concatenated_data_subset.h5ad', 'colon/data/04_integrated_data/dissected_mouse/'),
    'bone_marrow':   ('bone_marrow/data/03_concatenated_data/03_concatenated_data_subset.h5ad', 'bone_marrow/data/04_integrated_data/dissected_mouse/'),
    'fetal_liver':   ('fetal_liver/data/03_concatenated_data/03_concatenated_data_subset.h5ad', 'fetal_liver/data/04_integrated_data/dissected_mouse/'),
    'breast_cancer': ('breast_cancer/data/dissected_mouse/03_concatenated_data.h5ad', 'breast_cancer/data/dissected_mouse/'),
}
species_key = 'orig.ident'
out = 'figures/benchmark_summary/pca_variance.csv'

rows = []
for tissue, (f, d) in cases.items():
    a = sc.read_h5ad(f)
    cells = [c.strip() for c in open(d + 'eval_cells.txt')]
    a = a[a.obs_names.isin(cells)].copy()
    assert a.n_obs == len(cells), (tissue, a.n_obs, len(cells))
    layer = next(l for l in a.layers if l.startswith('data'))
    for sp in sorted(a.obs[species_key].astype(str).unique()):
        s = a[a.obs[species_key].astype(str) == sp].copy()
        s.X = s.layers[layer].copy()
        sc.pp.scale(s, max_value=10)
        X = np.asarray(s.X, dtype=np.float64)
        X = X[:, X.std(axis=0) > 0]                      # genes not expressed in this species carry no variance
        p = PCA(svd_solver='full').fit(X)
        r = p.explained_variance_ratio_
        rows.append(pd.DataFrame({'tissue': tissue, 'species': sp, 'pc': np.arange(1, len(r) + 1), 'var_ratio': r,
                                  'cum_ratio': np.cumsum(r), 'n_cells': s.n_obs, 'n_genes': X.shape[1]}))
        print(f'{tissue:14s} {sp:10s} {s.n_obs:7d} cells {X.shape[1]:4d} genes | PC1 {r[0]:.3f} | first 10 {r[:10].sum():.3f}, '
              f'20 {r[:20].sum():.3f}, 50 {r[:50].sum():.3f}', flush=True)
pd.concat(rows).to_csv(out, index=False)
print('done')
