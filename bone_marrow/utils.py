import os
import scanpy as sc
import pandas as pd
import scib
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score

BATCH = 'orig.ident'
# label set to score against. Default: EBCC for both species. A sensitivity run sets
# LABEL_SET=<name> and LABEL_COLS='hsapiens=<03 obs column>,mmusculus=<03 obs column>'
LABEL_SET = (os.environ.get('LABEL_SET') or 'ebcc')
LABEL_COLS = dict(x.split('=') for x in
                  (os.environ.get('LABEL_COLS') or 'hsapiens=ebcc_label,mmusculus=ebcc_label').split(','))
LABEL = 'ebcc_label' if LABEL_SET == 'ebcc' else f'label_{LABEL_SET}'

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

    # sensitivity label sets: build the label column per species from the 03 obs
    if LABEL not in raw.obs:
        species = raw.obs[BATCH].astype(str).values
        lab = pd.Series(index = raw.obs_names, dtype = object)
        for sp, col in LABEL_COLS.items():
            lab[species == sp] = raw.obs[col].astype(str).values[species == sp]
        raw.obs[LABEL] = lab.values

    # batch and labels always come from the unintegrated data
    for key in (BATCH, LABEL):
        adata_int.obs[key] = raw.obs[key].values
        adata_int.obs[key] = adata_int.obs[key].astype('category')

    # unintegrated reference for PCR: PCA on log-normalised data scaled across both species
    joint_pca(raw)

    is_embed = output_type == 'embed'

    if is_embed:
        sc.pp.neighbors(adata_int, use_rep = embed)
        adata_int.obsm["X_emb"] = adata_int.obsm[embed]

    # graph-only methods: pass the rep the method's own graph was built on. scib's
    # clustering rebuilds the kNN graph in place (from `embed`) whenever `embed` differs
    # from uns['neighbors']['params']['use_rep'], and every later metric then scores that
    # rebuilt graph instead of the method's (this silently replaced the SAMap, BBKNN UINMF
    # and BBKNN scVI-latent graphs with a plain kNN on X_pca)
    if not is_embed:
        embed = adata_int.uns.get('neighbors', {}).get('params', {}).get('use_rep', 'X_pca')
        # the graph LISI (C++) needs sorted CSR indices; align() reorders cells, which
        # can leave them unsorted ("malformed matrix format", no output, crash)
        for key in ('connectivities', 'distances'):
            if key in adata_int.obsp:
                g = adata_int.obsp[key].tocsr()
                g.sort_indices()
                adata_int.obsp[key] = g

    res = scib.me.metrics(
        raw, adata_int,
        batch_key = BATCH,
        label_key = LABEL,
        embed = embed,
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
    return res.rename(columns = {res.columns[0]: method})
