# Listing A.x: per-method metric computation (Python).
# Condensed from <tissue>/utils.py (score), sas_scores.py and extra_metrics.py;
# logging, file handling and per-tissue path code removed. Packages: scib, scanpy, scikit-learn.
import numpy as np, scanpy as sc, scib
from sklearn.neighbors import NearestNeighbors

BATCH, LABEL = 'orig.ident', 'ebcc_label'      # species; EBCC consensus cell-type label

# ---- 1. scib metrics -------------------------------------------------------------------------
def joint_pca(raw):
    """Unintegrated reference for PCR batch: log-normalised data, scaled across both species, PCA.
    The same PCA is the embedding of the 'Unintegrated' row that the improvement scores subtract."""
    raw.X = raw.layers['data'].copy()
    sc.pp.scale(raw, max_value=10)
    sc.pp.pca(raw)

def score(adata_int, raw, output_type, embed=None):
    """adata_int: one method's output; raw: unintegrated data, same cells in the same order.
    output_type: 'embed' (embedding methods) or 'knn' (graph-based: SAMap, BBKNN variants)."""
    for key in (BATCH, LABEL):                               # labels always from the unintegrated data
        adata_int.obs[key] = raw.obs[key].astype('category').values
    joint_pca(raw)
    is_embed = output_type == 'embed'
    if is_embed:
        sc.pp.neighbors(adata_int, use_rep=embed)           # 15-NN graph on the integrated embedding
        adata_int.obsm['X_emb'] = adata_int.obsm[embed]
    # graph-based methods are scored on their own saved graph (obsp['connectivities'])
    return scib.me.metrics(
        raw, adata_int, batch_key=BATCH, label_key=LABEL, embed=embed, type_=output_type,
        nmi_=True, ari_=True,                     # bio: NMI, ARI of Leiden clusters vs labels (scib defaults:
                                                  # Leiden at resolutions 0.2, 0.4, ..., 2.0; the clustering with
                                                  # the highest NMI is kept and ARI is computed on that same clustering;
                                                  # NMI with arithmetic averaging)
        graph_conn_=True,                         # batch: graph connectivity
        ilisi_=True, clisi_=True,                 # batch iLISI, bio cLISI (graph LISI, scaled 0-1)
        silhouette_=is_embed,                     # bio ASW_label, batch ASW_label/batch (embeddings only)
        pcr_=is_embed,                            # batch: PCR comparison with joint_pca(raw)
        isolated_labels_f1_=True, isolated_labels_asw_=is_embed,   # computed, not used in the score
        trajectory_=False, cell_cycle_=False, hvg_score_=False, kBET_=False)
# Breast cancer only: ASW_label and ASW_label/batch on a fixed 50,000-cell subsample
# (species proportions kept, seed 0), via scib.me.silhouette / scib.me.silhouette_batch.

# ---- 2. Mean average precision (MAP; bio) - scGLUE implementation ----------------------------
def mean_average_precision(x, y, neighbor_frac=0.01):
    """x: integrated embedding (cells x dims); y: cell-type labels."""
    k = max(round(len(y) * neighbor_frac), 1)                          # K = 1% of all cells
    nni = NearestNeighbors(n_neighbors=k + 1).fit(x).kneighbors(x, return_distance=False)
    match = y[nni[:, 1:]] == y[:, None]                                # drop self (first neighbour)
    prec = np.cumsum(match, axis=1) / np.arange(1, k + 1)              # precision at each rank
    n = match.sum(axis=1)
    ap = np.where(n > 0, (prec * match).sum(axis=1) / np.maximum(n, 1), 0.0)   # average precision
    return ap.mean()

# ---- 3. Seurat alignment score (SAS; batch) - scGLUE implementation --------------------------
def seurat_alignment_score(x, species, neighbor_frac=0.01, n_repeats=4, seed=20261030):
    rs = np.random.RandomState(seed)
    groups = [np.where(species == s)[0] for s in np.unique(species)]
    n, m = len(groups), min(g.size for g in groups)
    out = []
    for _ in range(n_repeats):                                         # subsample species to equal size
        sub = np.concatenate([rs.choice(g, m, replace=False) for g in groups])
        xs, ys = x[sub], species[sub]
        k = max(round(sub.size * neighbor_frac), 1)                    # K = 1% of the subsampled cells
        nni = NearestNeighbors(n_neighbors=k + 1).fit(xs).kneighbors(xs, return_distance=False)
        same = (ys[nni[:, 1:]] == ys[:, None]).sum(axis=1).mean()      # mean same-species neighbours
        out.append(min((k - same) * n / (k * (n - 1)), 1))           # = 1 - (same - k/n) / (k - k/n)
    return np.mean(out)
# MAP and SAS need an embedding: not defined for the graph-based methods (SAMap, BBKNN variants).

# ---- 4. LISI-CoM (reported only, not scored) -------------------------------------------------
def lisi_com(ilisi, clisi):
    return 2 * ilisi * clisi / (ilisi + clisi)     # harmonic mean of the scaled scib iLISI and cLISI
