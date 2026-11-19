# Values for the combined marker dot plot (analysis/marker_dotplot_figure.R), per tissue and species, by EBCC label:
#   python -I analysis/marker_dotplot_compute.py
# As Seurat's DotPlot (used for the per-dataset dot plots, source.R save_marker_dotplots): on the log-normalised layer,
# avg_exp = mean of expm1(value) over the cells of a label; pct = % of those cells with value > 0. (The per-gene scaling
# across labels is done in the figure script.) Genes: the markers on each panel, from figures/cross_tissue/marker_table.csv
# (analysis/marker_table.R, the lists used for the marker AUROC).
# Output: figures/cross_tissue/marker_dotplot_values.csv
# Optional: restrict to one half of the EBCC frac-selection split (set.seed(1); sample(n, floor(n / 2)) in R on the 02
# object = 'selection', the half f was tuned on; the rest = 'evaluation'). Drawn by calling R, so it is R's RNG:
#   python -I analysis/marker_dotplot_compute.py selection    -> marker_dotplot_values_selection.csv
# The selection half is checked against the EBCC label counts stored by the frac selection (<species>_ebcc_frac_auroc_per_class.csv).
import sys, os, re, subprocess
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp

# ---------------- config ----------------
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TISSUES = ['colon', 'bone_marrow', 'fetal_liver', 'breast_cancer']
LAYER = {('breast_cancer', 'human'): 'data', ('breast_cancer', 'mouse'): 'data_MERFISH'}   # default: data_Xenium
LABEL = 'ebcc_label'
MARKERS = f'{ROOT}/figures/cross_tissue/marker_table.csv'
HALF = sys.argv[1] if len(sys.argv) > 1 else 'all'      # all | selection | evaluation
OUT = f'{ROOT}/figures/cross_tissue/marker_dotplot_values' + ('' if HALF == 'all' else f'_{HALF}') + '.csv'
RSCRIPT = os.path.join(os.environ.get('CONDA_ENVS', os.path.expanduser('~/miniconda3/envs')), 'r45/bin/Rscript')
# -----------------------------------------
assert HALF in ('all', 'selection', 'evaluation'), HALF

def selection_mask(n):
    # R: set.seed(1); sample(n, floor(n / 2)) -> 1-based rows of the selection half
    env = dict(os.environ, R_LIBS_SITE = '', R_LIBS_USER = '')
    r = subprocess.run([RSCRIPT, '--vanilla', '-e', f'set.seed(1); cat(sample({n}, floor({n} / 2)), sep = "\\n")'],
                       capture_output = True, text = True, env = env, check = True)
    idx = np.array(r.stdout.split(), dtype = int) - 1
    m = np.zeros(n, bool); m[idx] = True
    assert m.sum() == n // 2
    return m

def check_selection(t, s, lab):
    # EBCC (chosen f) label counts on the selection half must equal those stored by the frac selection
    d = f'{ROOT}/{t}/data/02_annotated_data'
    f = pd.read_csv(f'{d}/{s}_ebcc_frac_auroc_chosen.csv')['frac'].iloc[0]
    pc = pd.read_csv(f'{d}/{s}_ebcc_frac_auroc_per_class.csv')
    pc = pc[pc.half.isin(['selection', 'select']) & np.isclose(pc.frac, f)]   # bone_marrow writes 'select'
    if pc.empty:
        sys.exit(f'{t} {s}: no stored selection-half counts for frac = {f}')
    key = lambda x: re.sub('[^a-z0-9]+', '_', str(x).lower())          # names compared as key() in analysis/vote_vs_ebcc.R
    got = pd.Series([key(l) for l in lab]).value_counts()
    bad = [(r.celltype, r.n_cells, got.get(key(r.celltype), 0)) for r in pc.itertuples() if got.get(key(r.celltype), 0) != r.n_cells]
    if bad:
        sys.exit(f'{t} {s}: selection half does not match the stored frac selection (type, stored, here): {bad}')
    return len(pc)

mt = pd.read_csv(MARKERS)
rows = []
for t in TISSUES:
    for s in ['human', 'mouse']:
        m = mt[(mt.tissue == t) & (mt.species == s) & mt.on_panel.notna()]
        genes = [(r.label, g.strip()) for r in m.itertuples() for g in r.on_panel.split(',')]
        a = ad.read_h5ad(f'{ROOT}/{t}/data/02_annotated_data/02_{s}_annotated.h5ad')
        layer = LAYER.get((t, s), 'data_Xenium')
        miss = sorted({g for _, g in genes if g not in a.var_names})
        if miss:
            sys.exit(f'{t} {s}: markers not in var_names: {miss}')
        uniq = list(dict.fromkeys(g for _, g in genes))
        X = a[:, uniq].layers[layer]
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        lab = a.obs[LABEL].astype(str).values
        if HALF != 'all':
            keep = selection_mask(a.n_obs)
            if HALF == 'evaluation': keep = ~keep
            else: n_checked = check_selection(t, s, lab[keep])
            X, lab = X[keep], lab[keep]
        for l in sorted(set(lab)):
            sub = X[lab == l]
            e = sub.copy(); e.data = np.expm1(e.data)
            avg = np.asarray(e.mean(axis = 0)).ravel()
            pct = 100 * np.asarray((sub > 0).mean(axis = 0)).ravel()
            rows.append(pd.DataFrame({'tissue': t, 'species': s, 'label': l, 'n_cells': sub.shape[0], 'gene': uniq,
                                      'avg_exp': avg, 'pct': pct}))
        # marker group of each gene (first group wins, as save_marker_dotplots)
        print(f'{t:14s} {s}: {len(lab)} of {a.n_obs} cells ({HALF}), {len(set(lab))} labels, {len(uniq)} marker genes (layer {layer})'
              + (f'; label counts match the stored selection half ({n_checked} types)' if HALF == 'selection' else ''))
out = pd.concat(rows)
gmap = pd.DataFrame([(r.tissue, r.species, r.label, g.strip(), i) for i, r in enumerate(mt[mt.on_panel.notna()].itertuples())
                     for g in r.on_panel.split(',')], columns = ['tissue', 'species', 'marker_of', 'gene', 'order'])
gmap = gmap.sort_values('order').drop_duplicates(['tissue', 'species', 'gene'])
out = out.merge(gmap.drop(columns = 'order'), on = ['tissue', 'species', 'gene'], how = 'left')
out.to_csv(OUT, index = False)
print(f'wrote {OUT} ({len(out)} rows)')
