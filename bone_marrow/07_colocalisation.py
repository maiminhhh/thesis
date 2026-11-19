print('1. Colocalisation - Setup')
# Spatial colocalisation of cell types in each species and its conservation across species.
#   human: its own EBCC labels (all cells)
#   mouse (dissected): labels transferred from human through the integration (06)
#   MOUSE_LABEL=ebcc_label python 07_colocalisation.py   # validation with the hidden mouse labels
# For each pair (a, b): among the k nearest neighbours of cells of type a, how often is the
# neighbour of type b, compared with label shuffles (positions fixed, type sizes fixed).
#   log2 ratio = log2(observed / mean shuffled)  -> effect size, comparable across species
#   z          = (observed - mean shuffled) / sd -> significance (grows with cell numbers)
# Conservation: correlation of the log2 ratios between species, overall and per cell type.
import os
N_CPU = 8
os.environ.setdefault("OMP_NUM_THREADS", str(N_CPU))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(N_CPU))
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

KS = [10, 30]          # neighbourhood sizes
N_PERM = 500
MIN_CELLS = 100        # a type is analysed only with at least this many cells in both species
SEED = 1
MOUSE_LABEL = os.environ.get('MOUSE_LABEL', 'celltype_transfer')
in_dir = 'data/06_label_transfer/dissected/'
out_dir = f'data/07_colocalisation/mouse_{MOUSE_LABEL}/'
fig_dir = f'../figures/bone_marrow/colocalisation/mouse_{MOUSE_LABEL}/'
os.makedirs(out_dir, exist_ok = True); os.makedirs(fig_dir, exist_ok = True)

cells = {'human': pd.read_csv(f'{in_dir}06_human_transfer_cells.csv').assign(label = lambda d: d.ebcc_label),
         'mouse': pd.read_csv(f'{in_dir}06_mouse_dissected_transfer_cells.csv').assign(
             label = lambda d: d[MOUSE_LABEL])}
for sp, d in cells.items():
    d['label'] = d.label.fillna('Unassigned').astype(str)
    print(f'{sp}: {len(d)} cells;', d.label.value_counts().to_dict())

counts = pd.DataFrame({sp: d.label.value_counts() for sp, d in cells.items()}).fillna(0)
types = sorted(counts.index[(counts.min(axis = 1) >= MIN_CELLS) & (counts.index != 'Unassigned')])
print(f'types with >= {MIN_CELLS} cells in both species: {types}')
counts.astype(int).to_csv(f'{out_dir}07_cell_counts.csv')

def colocalisation(d, k, rng):
    """observed neighbour fractions, shuffled mean / sd for every (a, b) pair of `types`"""
    xy = d[['x', 'y']].to_numpy()
    _, idx = NearestNeighbors(n_neighbors = k + 1, n_jobs = N_CPU).fit(xy).kneighbors(xy)
    idx = idx[:, 1:]                                            # drop the cell itself
    n = len(d)
    A = sparse.csr_matrix((np.ones(n * k), (np.repeat(np.arange(n), k), idx.ravel())), shape = (n, n))
    lab = pd.Categorical(d.label, categories = types)           # other labels -> NaN code (-1)
    code = lab.codes

    def fractions(code):
        L = sparse.csr_matrix((np.ones((code >= 0).sum()), (np.where(code >= 0)[0], code[code >= 0])),
                              shape = (n, len(types)))
        C = (L.T @ (A @ L)).toarray()                           # C[a, b] = neighbours of type b around type a
        return C / (np.asarray(L.sum(axis = 0)).ravel()[:, None] * k)
    obs = fractions(code)
    perm = np.stack([fractions(rng.permutation(code)) for _ in range(N_PERM)])
    return obs, perm.mean(0), perm.std(0)

def heatmaps(mats, title, out, vmax = 2, cbar = 'log2(observed / expected)', z = None):
    fig, axes = plt.subplots(1, len(mats), figsize = (5.2 * len(mats) + 1, 5))
    for ax, (name, m) in zip(np.atleast_1d(axes), mats.items()):
        im = ax.imshow(m.values, cmap = 'RdBu_r', vmin = -vmax, vmax = vmax)
        ax.set_xticks(range(len(types)), types, rotation = 45, ha = 'right', fontsize = 8)
        ax.set_yticks(range(len(types)), types, fontsize = 8)
        if z is not None:                                       # mark |z| > 3
            for i in range(len(types)):
                for j in range(len(types)):
                    if abs(z[name].values[i, j]) > 3:
                        ax.text(j, i, '*', ha = 'center', va = 'center', fontsize = 8)
        ax.set_title(name, fontsize = 10)
        ax.set_xlabel('neighbour type')
    np.atleast_1d(axes)[0].set_ylabel('centre type')
    fig.colorbar(im, ax = axes, label = cbar, fraction = 0.03)
    fig.suptitle(title, fontsize = 10)
    fig.savefig(out, dpi = 200, bbox_inches = 'tight')
    plt.close(fig)

summary = []
for k in KS:
    print(f'2. Colocalisation - k = {k}')
    ratio, z = {}, {}
    for sp, d in cells.items():
        obs, mu, sd = colocalisation(d, k, np.random.default_rng(SEED))
        ratio[sp] = pd.DataFrame(np.log2((obs + 1e-6) / (mu + 1e-6)), index = types, columns = types)
        z[sp] = pd.DataFrame((obs - mu) / np.where(sd > 0, sd, np.nan), index = types, columns = types)
        ratio[sp].to_csv(f'{out_dir}07_k{k}_{sp}_log2ratio.csv')
        z[sp].to_csv(f'{out_dir}07_k{k}_{sp}_z.csv')

    print(f'3. Colocalisation - conservation - k = {k}')
    h, m = ratio['human'].values.ravel(), ratio['mouse'].values.ravel()
    sig = (np.abs(z['human'].values.ravel()) > 3) | (np.abs(z['mouse'].values.ravel()) > 3)
    # same-type clustering (diagonal) is strong in most tissues and can carry the overall
    # correlation on its own, so between-type pairs are also scored separately
    off = ~np.eye(len(types), dtype = bool).ravel()
    summary.append(dict(k = k, n_types = len(types), pearson = pearsonr(h, m)[0], spearman = spearmanr(h, m)[0],
                        pearson_between = pearsonr(h[off], m[off])[0], spearman_between = spearmanr(h[off], m[off])[0],
                        pearson_homotypic = pearsonr(h[~off], m[~off])[0],
                        sign_agreement_sig_pairs = np.mean(np.sign(h[sig]) == np.sign(m[sig])),
                        n_sig_pairs = int(sig.sum())))
    # per cell type: how well its neighbourhood profile (row) is preserved
    per_type = pd.DataFrame({
        'pearson_row': [pearsonr(ratio['human'].loc[t], ratio['mouse'].loc[t])[0] for t in types],
        'pearson_row_between': [pearsonr(ratio['human'].loc[t].drop(t), ratio['mouse'].loc[t].drop(t))[0]
                                for t in types],
        'homotypic_human': np.diag(ratio['human']), 'homotypic_mouse': np.diag(ratio['mouse']),
        'n_human': counts.loc[types, 'human'].astype(int).values,
        'n_mouse': counts.loc[types, 'mouse'].astype(int).values}, index = types).round(3)
    per_type.to_csv(f'{out_dir}07_k{k}_conservation_per_type.csv')
    print(per_type.sort_values('pearson_row', ascending = False).to_string())

    heatmaps({'human (own EBCC labels)': ratio['human'], f'mouse ({MOUSE_LABEL})': ratio['mouse']},
             f'bone_marrow dissected: cell-type colocalisation, {k} nearest neighbours (* |z| > 3, {N_PERM} shuffles)',
             f'{fig_dir}heatmap_k{k}.png',
             z = {'human (own EBCC labels)': z['human'], f'mouse ({MOUSE_LABEL})': z['mouse']})

    fig, ax = plt.subplots(figsize = (5.5, 5))
    ax.scatter(h[off], m[off], s = 14, c = np.where(sig[off], 'black', 'lightgrey'), label = 'between types')
    ax.scatter(h[~off], m[~off], s = 30, marker = 'D', c = 'tab:red', label = 'same type (homotypic)')
    lim = np.nanmax(np.abs(np.r_[h, m])) * 1.05
    ax.plot([-lim, lim], [-lim, lim], ls = '--', c = 'grey', lw = 0.8)
    ax.axhline(0, c = 'grey', lw = 0.5); ax.axvline(0, c = 'grey', lw = 0.5)
    ax.set_xlabel('human: log2(observed / expected)'); ax.set_ylabel('mouse: log2(observed / expected)')
    ax.set_title(f'k = {k}: Pearson r = {summary[-1]["pearson"]:.2f} all pairs, '
                 f'{summary[-1]["pearson_between"]:.2f} between types\n({len(types)} types; grey = |z| <= 3 in both)',
                 fontsize = 9)
    ax.legend(fontsize = 7, loc = 'lower right')
    fig.savefig(f'{fig_dir}scatter_k{k}.png', dpi = 200, bbox_inches = 'tight')
    plt.close(fig)

pd.DataFrame(summary).round(3).to_csv(f'{out_dir}07_conservation_summary.csv', index = False)
print(pd.DataFrame(summary).round(3).to_string(index = False))
print('Colocalisation - Done!')
