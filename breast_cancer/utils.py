import os
import anndata as ad
import numpy as np
import scanpy as sc
import pandas as pd
import scib
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score

BATCH = 'orig.ident'
LABEL = 'ebcc_label'

# label sensitivity check (../EBCC_PIPELINE.md section 9b). Unset = the primary EBCC labels, unchanged.
# Otherwise LABEL_SET names the run and LABEL_COLS='hsapiens=<col>,mmusculus=<col>' gives, per species,
# the 03 obs column the labels are taken from.
LABEL_SET = os.environ.get('LABEL_SET', 'ebcc')
LABEL_COLS = None
if LABEL_SET != 'ebcc':
    if 'LABEL_COLS' not in os.environ:
        raise SystemExit(f'LABEL_SET={LABEL_SET} needs LABEL_COLS=hsapiens=<col>,mmusculus=<col>')
    LABEL_COLS = dict(kv.split('=') for kv in os.environ['LABEL_COLS'].split(','))

# the label of every cell under LABEL_COLS, taken per species from the 03 data
def labels(raw):
    species = raw.obs[BATCH].astype(str).values
    lab = pd.Series(index=raw.obs_names, dtype=object)
    for sp, col in LABEL_COLS.items():
        lab[species == sp] = raw.obs.loc[species == sp, col].astype(str).values
    if lab.isna().any():
        raise ValueError(f'{lab.isna().sum()} cells have no label: LABEL_COLS must cover {sorted(set(species))}')
    return lab.values

# SAMap keeps the species objects' cell names; add the species prefix only if it is missing
def samap_obs_names(adata):
    species = ['mmusculus' if s == 'mm' else 'hsapiens' for s in adata.obs['orig.ident'].astype(str)]
    return pd.Index([c if c.startswith(f'{s}_') else f'{s}_{c}'
                     for s, c in zip(species, adata.obs_names)])

# evaluation cells written by 04.11: cells of the concatenated data whose UINMF H row has
# >= 2 non-zero factors. Every method and both baselines are scored on exactly these cells.
def eval_cells(g):
    path = f'data/{g}/eval_cells.txt'
    if not os.path.exists(path):
        raise FileNotFoundError(f'{path} is missing: run 04.11 (LIGER UINMF) for {g} first')
    with open(path) as f:
        return pd.Index([l.strip() for l in f if l.strip()])

# keep the evaluation cells of dataset g that the integrated output has, in raw's order
def align(adata_int, raw, g):
    cells = eval_cells(g)
    common = raw.obs_names[raw.obs_names.isin(cells) & raw.obs_names.isin(adata_int.obs_names)]
    print(f'{g}: scoring {len(common)} cells | evaluation cells: {len(cells)} | '
          f'missing from the integrated output: {len(cells) - len(common)}')
    return adata_int[common].copy(), raw[common].copy()

# baseline PCA for PCR_batch (../EBCC_PIPELINE.md section 7): the 03 scale.data layer is scaled within each
# species, which removes most of the species effect, so scale the log-normalised data across both species here
def joint_pca(adata):
    data_layer = next((l for l in adata.layers if l.startswith('data')), None)
    if data_layer is not None:
        adata.X = adata.layers[data_layer].copy()
    sc.pp.scale(adata, max_value = 10)
    sc.pp.pca(adata)

# scoring function for scib
def score(adata_int, method, output_type, raw, embed = None):
    """
    Score 1 integration result with scib
    output_type: 'embed' for embedding methods, 'knn' for graph-only methods
    raw: the unintegrated data (03_concatenated_data.h5ad) for the same cells
    """

    assert (adata_int.obs_names == raw.obs_names).all(), 'cells must match and be in the same order'

    # batch and labels always come from the unintegrated data
    if LABEL_COLS is not None:
        raw.obs[LABEL] = pd.Categorical(labels(raw))
    for key in (BATCH, LABEL):
        adata_int.obs[key] = raw.obs[key].values
        adata_int.obs[key] = adata_int.obs[key].astype('category')

    # unintegrated reference: log-normalised data, PCA on the genes actually used
    # (X is raw counts, and the stored X_pca was computed before genes were dropped)
    joint_pca(raw)

    is_embed = output_type == 'embed'

    if is_embed:
        sc.pp.neighbors(adata_int, use_rep = embed)
        adata_int.obsm["X_emb"] = adata_int.obsm[embed]
    else:
        # align() reorders cells, which leaves the stored graph's column indices unsorted; scib's LISI binary
        # rejects unsorted input without raising, so sort them (values are unchanged)
        for k in adata_int.obsp:
            adata_int.obsp[k].sort_indices()

    res = scib.me.metrics(
        raw, adata_int,
        batch_key = BATCH,
        label_key = LABEL,
        embed = embed if is_embed else None,   # knn: score the stored graph
        type_ = output_type,

        # graph-based metrics: all methods
        nmi_ = True,
        ari_ = True,
        isolated_labels_f1_ = True,
        graph_conn_ = True,
        ilisi_ = True, clisi_ = True,

        # embedding-based metrics: embedding methods only
        # (the silhouette metrics are computed below on a subsample)
        silhouette_ = False,
        isolated_labels_asw_ = False,
        pcr_ = is_embed,

        # switched off
        trajectory_ = False,
        cell_cycle_ = False,
        hvg_score_ = False,
        kBET_ = False       # calculate in R later
    )

    # silhouettes compare every pair of cells, so on large datasets they use a fixed subsample
    if is_embed:
        sub = adata_int[asw_cells(adata_int.obs_names, adata_int.obs[BATCH])].copy()
        print(f'ASW metrics on {sub.n_obs} / {adata_int.n_obs} cells')
        res.loc['ASW_label'] = scib.me.silhouette(sub, LABEL, embed)
        res.loc['ASW_label/batch'] = scib.me.silhouette_batch(sub, BATCH, LABEL, embed, verbose = False)
        res.loc['isolated_label_silhouette'] = scib.me.isolated_labels_asw(sub, LABEL, BATCH, embed, verbose = False)
    return res.rename(columns = {res.columns[0]: method})

# fixed subsample for the silhouette metrics: ASW_N cells, keeping the species proportions;
# drawn from the sorted cell names, so every method on the same cells gets the same subsample
ASW_N = 50000

def asw_cells(obs_names, batch, n = None, seed = 0):
    n = ASW_N if n is None else n
    if len(obs_names) <= n:
        return obs_names
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({'cell': obs_names, 'batch': pd.Series(batch).astype(str).values}).sort_values('cell')
    picked = []
    for _, d in df.groupby('batch', sort = True):
        picked += list(rng.choice(d['cell'].to_numpy(), round(n * len(d) / len(df)), replace = False))
    return obs_names[obs_names.isin(picked)]
