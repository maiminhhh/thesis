print('1. Label transfer on scVI - Setup')
# kNN label transfer across species in the scVI latent space, then compare the
# transferred labels with each cell's own (EBCC) label.
#   python 06_label_transfer_scVI.py                  # human -> mouse, dissected
#   python 06_label_transfer_scVI.py genes2           # other datasets
#   DIRECTION=mouse_to_human python 06_label_transfer_scVI.py
# Default human -> mouse: the human labels cover every cell type (the dissected mouse
# has 17 T cells and 7 plasma cells, so a mouse reference can barely call those).
import os
N_CPU = 8
os.environ.setdefault("OMP_NUM_THREADS", str(N_CPU))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(N_CPU))
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from datasets import DATASETS

K = 30
LABEL = 'ebcc_label'
DIRECTION = os.environ.get('DIRECTION', 'human_to_mouse')
REF, QUERY = {'human_to_mouse': ('hsapiens', 'mmusculus'),
              'mouse_to_human': ('mmusculus', 'hsapiens')}[DIRECTION]
NAMES = {'hsapiens': 'human', 'mmusculus': 'mouse'}
names = sys.argv[1:] or ['dissected']

def heatmap(tab, title, out, fmt = '{:.0f}', cbar = '% of row'):
    fig, ax = plt.subplots(figsize = (1 + 0.6 * tab.shape[1], 1.2 + 0.5 * tab.shape[0]))
    im = ax.imshow(tab.values, cmap = 'Blues', vmin = 0, vmax = 100, aspect = 'auto')
    ax.set_xticks(range(tab.shape[1]), tab.columns, rotation = 45, ha = 'right')
    ax.set_yticks(range(tab.shape[0]), tab.index)
    for i in range(tab.shape[0]):
        for j in range(tab.shape[1]):
            v = tab.values[i, j]
            if v >= 0.5:
                ax.text(j, i, fmt.format(v), ha = 'center', va = 'center', fontsize = 7,
                        color = 'white' if v > 60 else 'black')
    ax.set_xlabel(f'after: label transferred from {NAMES[REF]} (scVI, {K}-NN)')
    ax.set_ylabel(f'before: {NAMES[QUERY]} {LABEL}')
    ax.set_title(title, fontsize = 10)
    fig.colorbar(im, ax = ax, label = cbar, fraction = 0.04)
    fig.tight_layout()
    fig.savefig(out, dpi = 200)
    plt.close(fig)

for name in names:
    d = DATASETS[name]
    print(f'2. Label transfer {NAMES[REF]} -> {NAMES[QUERY]} - {name}')
    adata = sc.read_h5ad(f"{d['out']}04.07_scVI.h5ad")
    sp = adata.obs['orig.ident'].astype(str).values
    lab = adata.obs[LABEL].astype(str).values
    Z = adata.obsm['X_scvi']
    ref, query = sp == REF, sp == QUERY

    nn = NearestNeighbors(n_neighbors = K, n_jobs = N_CPU).fit(Z[ref])
    _, idx = nn.kneighbors(Z[query])
    votes = pd.DataFrame(lab[ref][idx])
    transferred = votes.mode(axis = 1)[0].values
    confidence = (votes.values == transferred[:, None]).mean(axis = 1)

    res = pd.DataFrame({'before': lab[query], 'after': transferred, 'confidence': confidence},
                       index = adata.obs_names[query])
    out_dir = f'data/06_label_transfer/{name}/'
    fig_dir = f'../figures/bone_marrow/label_transfer/'
    os.makedirs(out_dir, exist_ok = True); os.makedirs(fig_dir, exist_ok = True)
    tag = f'{DIRECTION}_{name}'
    res.to_csv(f'{out_dir}06_scVI_{DIRECTION}_labels.csv')
    pd.Series(adata.obs_names[ref]).to_csv(f'{out_dir}06_scVI_{DIRECTION}_reference_cells.txt',
                                           index = False, header = False)

    print(f'3. Label transfer - Compare before / after - {name}')
    order = sorted(set(lab))
    counts = pd.crosstab(res.before, res.after).reindex(index = order, columns = order, fill_value = 0)
    counts = counts.loc[counts.sum(axis = 1) > 0]
    pct = 100 * counts.div(counts.sum(axis = 1), axis = 0)
    counts.to_csv(f'{out_dir}06_scVI_{DIRECTION}_confusion_counts.csv')
    agree = (res.before == res.after).mean()
    per_type = pd.DataFrame({'n': counts.sum(axis = 1),
                             'kept_pct': np.diag(pct.reindex(columns = pct.index)).round(1),
                             'top_other': [pct.loc[t].drop(t).idxmax() for t in pct.index],
                             'top_other_pct': [pct.loc[t].drop(t).max().round(1) for t in pct.index],
                             'median_confidence': res.groupby('before').confidence.median().reindex(pct.index).round(2)})
    per_type.to_csv(f'{out_dir}06_scVI_{DIRECTION}_per_type.csv')
    print(f'{name}: {query.sum()} {NAMES[QUERY]} cells, overall agreement {100 * agree:.1f}%, '
          f'median confidence {np.median(confidence):.2f}')
    print(per_type.to_string())

    # composition before vs after (does the transfer reshape the cell-type mix?)
    comp = pd.DataFrame({'before': res.before.value_counts(normalize = True),
                         'after': res.after.value_counts(normalize = True)}).reindex(order).fillna(0)
    (100 * comp).round(2).to_csv(f'{out_dir}06_scVI_{DIRECTION}_composition.csv')

    heatmap(pct, f'bone_marrow {name}: {NAMES[QUERY]} labels before vs after transfer from '
                 f'{NAMES[REF]} on scVI\n(row % of {NAMES[QUERY]} cells; overall agreement {100 * agree:.1f}%)',
            f'{fig_dir}heatmap_{tag}.png')
    print(f'Label transfer - {name} - Done!')
