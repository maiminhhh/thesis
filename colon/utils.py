import os
import scanpy as sc
import pandas as pd
import scib
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score

BATCH = 'orig.ident'

# label set the bio metrics use (section 9b of ../EBCC_PIPELINE.md). EBCC is the primary label set:
# LABEL_SET unset (or 'ebcc') scores on ebcc_label exactly as before. A robustness rescoring sets
# LABEL_SET=<name> and LABEL_COLS='hsapiens=<col>,mmusculus=<col>' (each species' 03 obs column)
LABEL_SET = os.environ.get('LABEL_SET') or 'ebcc'
if LABEL_SET == 'ebcc':
    LABEL, LABEL_COLS = 'ebcc_label', None
else:
    LABEL = f'label_{LABEL_SET}'
    LABEL_COLS = dict(kv.strip().split('=', 1) for kv in os.environ['LABEL_COLS'].split(','))
    assert set(LABEL_COLS) == {'hsapiens', 'mmusculus'}, f'LABEL_COLS needs hsapiens= and mmusculus=: {LABEL_COLS}'

# the LABEL column for a non-default label set: each species' cells take their own column.
# 'colA|colB' falls back to colB where colA is missing (SingleR's pruned calls -> same-species Seurat)
def set_labels(obs):
    lab = pd.Series(pd.NA, index = obs.index, dtype = object)
    for species, spec in LABEL_COLS.items():
        is_sp = (obs[BATCH] == species).to_numpy()
        for col in spec.split('|'):
            fill = is_sp & lab.isna().to_numpy()
            lab[fill] = obs.loc[fill, col].astype(object).to_numpy()
    assert lab.notna().all(), f'{LABEL_SET}: {lab.isna().sum()} cells without a label'
    obs[LABEL] = lab.astype('category')

# the cells every method is scored on, written by 04.11 (cells with a usable UINMF H row)
def eval_cells(out_dir):
    path = f'{out_dir}eval_cells.txt'
    if not os.path.exists(path):
        raise FileNotFoundError(f'{path} not found: run 04.11_integrate_LIGER_UINMF.R first')
    return pd.Index(pd.read_csv(path, header=None)[0].astype(str))

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

# keep only the cells present in both objects (and in `cells`, if given), in raw's order
def align(adata_int, raw, cells = None):
    common = raw.obs_names.intersection(adata_int.obs_names, sort=False)
    if cells is not None:
        common = common.intersection(cells, sort=False)
    print(f'common: {len(common)} | only in integrated: {adata_int.n_obs - len(common)} '
          f'| only in raw: {raw.n_obs - len(common)}')
    if cells is not None and len(common) < len(cells):
        print(f'warning: {len(cells) - len(common)} evaluation cells missing from this method')
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
    if LABEL_SET != 'ebcc':
        set_labels(raw.obs)
    for key in (BATCH, LABEL):
        adata_int.obs[key] = raw.obs[key].values
        adata_int.obs[key] = adata_int.obs[key].astype('category')

    # unintegrated reference for PCR: PCA on log-normalised data scaled across both species
    joint_pca(raw)

    is_embed = output_type == 'embed'

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
        # SAMap's uns['neighbors'] holds only 'params': fill in whatever keys are missing
        nb = dict(adata_int.uns.get('neighbors', {}))
        nb.setdefault('connectivities_key', 'connectivities')
        nb.setdefault('distances_key', 'distances')
        nb.setdefault('params', {'method': 'umap'})
        adata_int.uns['neighbors'] = nb

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
        kBET_ = False       # calculate in R later
    )
    return res.rename(columns = {res.columns[0]: method})
