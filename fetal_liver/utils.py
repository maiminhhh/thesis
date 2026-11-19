import scanpy as sc
import pandas as pd
import scib
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score

import os
BATCH = 'orig.ident'
LABEL = 'ebcc_label'
# labels for the robustness check (../EBCC_PIPELINE.md section 9b): with LABEL_SET set,
# LABEL_COLS='hsapiens=<col>,mmusculus=<col>' gives each species' 03 obs column.
# Both unset = the EBCC labels, unchanged
LABEL_SET = os.environ.get('LABEL_SET') or 'ebcc'
LABEL_COLS = os.environ.get('LABEL_COLS') or ''
assert (LABEL_SET == 'ebcc') == (LABEL_COLS == ''), 'set LABEL_SET and LABEL_COLS together'

# with LABEL_COLS: put each species' column into raw.obs[LABEL]
def set_labels(raw):
    if LABEL_SET == 'ebcc':
        return
    cols = dict(kv.split('=') for kv in LABEL_COLS.split(','))
    species = raw.obs[BATCH].astype(str)
    assert set(species) <= set(cols), f'LABEL_COLS {cols} misses species {set(species) - set(cols)}'
    lab = pd.Series(index = raw.obs_names, dtype = object)
    for sp, col in cols.items():
        lab[species == sp] = raw.obs.loc[species == sp, col].astype(str)
    assert lab.notna().all(), 'cells without a label'
    raw.obs[LABEL] = lab.values
    print(f'labels: {LABEL_SET} {cols}')

# baseline PCA: log-normalised data scaled jointly across both species, then PCA.
# The 03 scale.data layer is scaled within each species, which removes most of the
# species effect (PCR_batch baseline ~0, so most methods scored 0).
# X is raw counts, and the stored X_pca was computed before genes were dropped.
def joint_pca(adata):
    data_layer = next((l for l in adata.layers if l.startswith('data')), None)
    if data_layer is not None:
        adata.X = adata.layers[data_layer].copy()
    sc.pp.scale(adata, max_value = 10)
    sc.pp.pca(adata)

# the cells every method and baseline is scored on (written by 04.11): in the
# concatenated data, with >= 2 non-zero UINMF factors
def eval_cells(out):
    with open(f'{out}eval_cells.txt') as f:
        return pd.Index([line.strip() for line in f if line.strip()])

# keep only the evaluation cells, in raw's order; every one must be in both objects
def align(adata_int, raw, cells):
    missing = cells.difference(adata_int.obs_names).union(cells.difference(raw.obs_names))
    assert len(missing) == 0, f'{len(missing)} evaluation cells missing, e.g. {list(missing[:3])}'
    common = raw.obs_names.intersection(cells, sort=False)
    print(f'evaluation cells: {len(common)} | integrated has {adata_int.n_obs}, raw has {raw.n_obs}')
    return adata_int[common].copy(), raw[common].copy()

# scoring function for scib
def score(adata_int, method, output_type, raw, embed = None):
    """
    Score 1 integration result with scib
    output_type: 'embed' for embedding methods, 'knn' for graph-only methods
    raw: the unintegrated data (03_concatenated_data.h5ad) for the same cells
    """

    assert (adata_int.obs_names == raw.obs_names).all(), 'cells must match and be in the same order'

    # batch and labels always come from the unintegrated data
    set_labels(raw)
    for key in (BATCH, LABEL):
        adata_int.obs[key] = raw.obs[key].values
        adata_int.obs[key] = adata_int.obs[key].astype('category')

    is_embed = output_type == 'embed'

    # diagnostic only: PCR_batch against the old baseline (03 scale.data, scaled within each
    # species), kept to report how the fix changed the scores; not used in any score
    if is_embed:
        old = raw.copy()
        data_layer = next((l for l in old.layers if l.startswith('data')), None)
        scale_layer = next((l for l in old.layers if l.startswith('scale.data')), None)
        if data_layer is not None:
            old.X = old.layers[data_layer]
        sc.pp.pca(old, layer = scale_layer)
        pcr_within_species = scib.me.pcr_comparison(old, adata_int, covariate = BATCH, embed = embed)
        del old

    # unintegrated reference for PCR_batch: PCA of the data scaled across both species
    joint_pca(raw)

    if is_embed:
        sc.pp.neighbors(adata_int, use_rep = embed)
        adata_int.obsm["X_emb"] = adata_int.obsm[embed]
    else:
        # graph-only methods (BBKNN, SAMap) are scored on their own kNN graph: passing an
        # embedding here makes scib rebuild the neighbour graph from it, discarding the
        # integrated graph. Register the graph if the method did not (e.g. SAMap)
        assert 'connectivities' in adata_int.obsp, 'graph method without obsp["connectivities"]'
        # scib's graph-LISI C++ reader needs CSR with sorted column indices (BBKNN's are not)
        conn = adata_int.obsp['connectivities'].tocsr()
        conn.sum_duplicates()
        conn.sort_indices()
        adata_int.obsp['connectivities'] = conn
        if 'neighbors' not in adata_int.uns:
            adata_int.uns['neighbors'] = {'connectivities_key': 'connectivities',
                                          'distances_key': 'distances',
                                          'params': {'method': 'umap'}}

    res = scib.me.metrics(
        raw, adata_int,
        batch_key = BATCH,
        label_key = LABEL,
        embed = embed if is_embed else None,
        type_ = output_type,

        # graph-based metrics: all methods
        nmi_ = True,
        ari_ = True,
        isolated_labels_f1_ = True,
        graph_conn_ = True,
        ilisi_ = True, clisi_ = True,

        # embedding-based metrics: embedding methods only
        silhouette_ = is_embed,
        isolated_labels_asw_ = is_embed,
        pcr_ = is_embed,

        # switched off
        trajectory_ = False,
        cell_cycle_ = False,
        hvg_score_ = False,
        kBET_ = False
    )
    if is_embed:
        res.loc['PCR_batch_within_species'] = pcr_within_species
    return res.rename(columns = {res.columns[0]: method})
