# Bone marrow (dissected) UMAPs with identical settings for the appendix method-comparison figure (analysis/bm_umap_figure.R).
# Nothing is integrated or refitted: each method's saved output is read and only the UMAP is (re)computed.
#   embedding methods (unintegrated joint PCA, scANVI, Harmony on scVI-latent): kNN graph from the embedding
#     (sc.pp.neighbors, settings below), then sc.tl.umap
#   graph-based methods (SAMap, BBKNN): sc.tl.umap directly on the method's own saved graph (obsp['connectivities'])
# UMAP is run on all cells of each file (so a graph is not cut); the figure shows the evaluation cells present in every method.
# Every UMAP starts from a random layout (init_pos='random'): spectral initialisation failed for three of the five graphs, so
# random is used for all to keep the settings identical. analysis/appendix_unintegrated_umap_compute.py imports these settings.
# Usage: python analysis/bm_umap_compute.py
#   -> figures/benchmark_summary/bm_umap/<method>.csv (cell, umap_1, umap_2), settings.txt
import os, h5py, numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp

# ---- config ----
base = 'bone_marrow/data/04_integrated_data/dissected_mouse/'
concat = 'bone_marrow/data/03_concatenated_data/03_concatenated_data_subset.h5ad'
out = 'figures/benchmark_summary/bm_umap'
methods = {  # name: (file, embedding key or None for a graph-based method)
    'Unintegrated':        (concat, 'X_pca'),
    'scANVI':              (base + '04.08_scANVI.h5ad', 'X_scanvi'),
    'SAMap':               (base + '04.02_SAMap.h5ad', None),
    'BBKNN':               (base + '04.09_BBKNN.h5ad', None),
    'Harmony scVI-latent': (base + '04.21_Harmony_scVI.h5ad', 'X_harmony'),
}
n_neighbors = 15          # kNN graph of embedding methods (scanpy default; = SAMap's saved graph)
metric = 'euclidean'
min_dist = 0.5            # sc.tl.umap defaults
spread = 1.0
seed = 0
init_pos = 'random'      # identical for every UMAP (spectral failed for SAMap, BBKNN, Harmony scVI-latent)

def cells_of(h):
    return h['obs'][h['obs'].attrs['_index']][:].astype(str)

def csr(g):
    return sp.csr_matrix((g['data'][:], g['indices'][:], g['indptr'][:]), shape=tuple(g.attrs['shape']))

def umap_of(cells, emb=None, graph=None, rep='X_emb'):
    """UMAP with the settings above, from an embedding (kNN graph built here) or a saved connectivity graph."""
    a = ad.AnnData(obs=pd.DataFrame(index=cells))
    if emb is not None:
        a.obsm[rep] = emb
    else:
        a.obsp['connectivities'] = graph
    if emb is not None:
        sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep=rep, metric=metric, random_state=seed)
    else:
        a.uns['neighbors'] = {'connectivities_key': 'connectivities', 'distances_key': 'distances',
                              'params': {'method': 'umap', 'metric': metric, 'n_neighbors': n_neighbors}}
    sc.tl.umap(a, min_dist=min_dist, spread=spread, random_state=seed, init_pos=init_pos)
    return a.obsm['X_umap']

def settings_text():
    return (f'scanpy {sc.__version__}; embedding methods: sc.pp.neighbors(n_neighbors={n_neighbors}, metric={metric}, '
            f'random_state={seed}, all embedding dimensions); graph-based methods: their saved connectivities; '
            f'sc.tl.umap(min_dist={min_dist}, spread={spread}, random_state={seed}, init_pos={init_pos}) for all.')

if __name__ == '__main__':
  os.makedirs(out, exist_ok=True)
  for name, (f, rep) in methods.items():
    with h5py.File(f, 'r') as h:
        cells = cells_of(h)
        u = umap_of(cells, emb=h['obsm/' + rep][:], rep=rep) if rep is not None else umap_of(cells, graph=csr(h['obsp/connectivities']))
    pd.DataFrame({'cell': cells, 'umap_1': u[:, 0], 'umap_2': u[:, 1]}).to_csv(
        os.path.join(out, name.replace(' ', '_') + '.csv'), index=False)
    print(name, len(cells), 'cells', 'embedding ' + rep if rep else 'own graph', flush=True)

  with open(os.path.join(out, 'settings.txt'), 'w') as fh:
      fh.write(settings_text() + '\n')
  print('done')
