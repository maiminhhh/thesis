# datasets every 04 integration / 05 evaluation script runs on; each is a folder data/<name>/
# holding 02_human_annotated, 02_mouse_annotated and 03_concatenated_data (.h5ad). Mirrors integration_sets.R.
import sys

DATASETS = ['full_mouse', 'dissected_mouse', 'genes2', 'genes4', 'genes6', 'genes8']


def selected():
    """datasets named on the command line, or all of them
    (e.g. `python 04.05_integrate_Scanorama.py genes2 genes4`)"""
    names = sys.argv[1:] or DATASETS
    unknown = [n for n in names if n not in DATASETS]
    if unknown:
        raise SystemExit(f'unknown datasets: {unknown}; choose from {DATASETS}')
    return names


def scores_dir(name):
    """folder for a dataset's 05 scores: data/<name>/ for the primary EBCC labels, and
    data/<name>/label_<set>/ for a LABEL_SET run of the label sensitivity check"""
    import os
    label_set = os.environ.get('LABEL_SET', 'ebcc')
    out = f'data/{name}/' if label_set == 'ebcc' else f'data/{name}/label_{label_set}/'
    os.makedirs(out, exist_ok=True)
    return out
