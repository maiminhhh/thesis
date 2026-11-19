# datasets every 04 script runs on: comment out an entry to skip it (mirrors integration_datasets.R)
#   concat = concatenated input from 03 (no extension: .rds for R, .h5ad for python)
#   inputs = folder with the per-species 02_{human,mouse}_annotated objects (Ensembl ids)
#   out    = output folder
#   scores = folder the 05 evaluation scores are written to
import os
import numpy as np
import pandas as pd

DATASETS = {
    'full_mouse':      dict(concat = 'data/03_concatenated_data/03_concatenated_data_full',
                            inputs = 'data/03_concatenated_data/full_mouse/',
                            out    = 'data/04_integrated_data/full_mouse/',
                            scores = 'data/05_evaluated_score/full_mouse/'),
    'dissected_mouse': dict(concat = 'data/03_concatenated_data/03_concatenated_data_subset',
                            inputs = 'data/03_concatenated_data/dissected_mouse/',
                            out    = 'data/04_integrated_data/dissected_mouse/',
                            scores = 'data/05_evaluated_score/dissected_mouse/'),
    'genes2':          dict(concat = 'data/genes2/03_concatenated_data',
                            inputs = 'data/genes2/',
                            out    = 'data/genes2/',
                            scores = 'data/genes2/'),
    'genes4':          dict(concat = 'data/genes4/03_concatenated_data',
                            inputs = 'data/genes4/',
                            out    = 'data/genes4/',
                            scores = 'data/genes4/'),
    'genes6':          dict(concat = 'data/genes6/03_concatenated_data',
                            inputs = 'data/genes6/',
                            out    = 'data/genes6/',
                            scores = 'data/genes6/'),
    'genes8':          dict(concat = 'data/genes8/03_concatenated_data',
                            inputs = 'data/genes8/',
                            out    = 'data/genes8/',
                            scores = 'data/genes8/'),
}
# 05 with LABEL_SET / LABEL_COLS (utils.py) scores into a label_<set>/ subfolder
_label_set = os.environ.get('LABEL_SET') or 'ebcc'
if _label_set != 'ebcc':
    for d in DATASETS.values():
        d['scores'] = f"{d['scores']}label_{_label_set}/"
for d in DATASETS.values():
    for p in (d['out'], d['scores']):
        os.makedirs(p, exist_ok = True)

# dataset names on the command line run only those, e.g. python 04.07_integrate_scVI.py genes2 genes4
import sys
_run_only = sys.argv[1:]
if _run_only:
    assert all(n in DATASETS for n in _run_only), f'unknown dataset in {_run_only}'
    DATASETS = {n: DATASETS[n] for n in _run_only}

# UINMF H rows for the cells in adata, dropping cells with no H row or an
# all-zero / non-finite one (same rule as the R helpers in source.R)
def load_uinmf_H(adata, out):
    H = pd.read_csv(f'{out}uinmf_H_raw.csv', index_col = 0)
    usable = np.isfinite(H.to_numpy()).all(axis = 1) & (H.abs().sum(axis = 1).to_numpy() > 0)
    common = adata.obs_names.intersection(H.index[usable], sort = False)
    print(f'UINMF: keeping {len(common)} / {adata.n_obs} cells '
          f'({adata.n_obs - len(common)} dropped: missing or all-zero H)')
    adata = adata[common].copy()
    adata.obsm['X_uinmf_H'] = H.loc[adata.obs_names].to_numpy(dtype = float)
    return adata
