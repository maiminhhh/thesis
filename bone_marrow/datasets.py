# datasets every 04 integration / 05 evaluation script runs on: concatenated input,
# the per-species inputs (SAMap, Ensembl IDs), the 04 output folder and the 05
# score folder. Mirrors integration_sets.R.
DATASETS = {
    'full': {
        'concat': 'data/03_concatenated_data/03_concatenated_data_full.h5ad',
        'human':  'data/03_concatenated_data/03_human.h5ad',
        'mouse':  'data/03_concatenated_data/03_mouse_full.h5ad',
        'out':    'data/04_integrated_data/full_mouse/',
        'scores': 'data/05_evaluated_score/full_mouse/'},
    'dissected': {
        'concat': 'data/03_concatenated_data/03_concatenated_data_subset.h5ad',
        'human':  'data/03_concatenated_data/03_human.h5ad',
        'mouse':  'data/03_concatenated_data/03_mouse_subset.h5ad',
        'out':    'data/04_integrated_data/dissected_mouse/',
        'scores': 'data/05_evaluated_score/dissected_mouse/'},
}
for g in [2, 4, 6, 8]:
    DATASETS[f'genes{g}'] = {
        'concat': f'data/genes{g}/03_concatenated_data.h5ad',
        'human':  f'data/genes{g}/02_human_annotated.h5ad',
        'mouse':  f'data/genes{g}/02_mouse_annotated.h5ad',
        'out':    f'data/genes{g}/',
        'scores': f'data/genes{g}/'}


# sensitivity runs (LABEL_SET other than ebcc) write scores to a label_<set>/ subfolder
import os as _os
_label_set = (_os.environ.get('LABEL_SET') or 'ebcc')
if _label_set != 'ebcc':
    for _d in DATASETS.values():
        _d['scores'] = f"{_d['scores']}label_{_label_set}/"
        _os.makedirs(_d['scores'], exist_ok = True)


def selected():
    """(name, paths) for the datasets named on the command line, or all of them
    (e.g. `python 04.05_integrate_Scanorama.py genes2 genes4`)"""
    import sys
    names = sys.argv[1:] or list(DATASETS)
    unknown = [n for n in names if n not in DATASETS]
    if unknown:
        raise SystemExit(f'unknown datasets: {unknown}; choose from {list(DATASETS)}')
    return [(n, DATASETS[n]) for n in names]
