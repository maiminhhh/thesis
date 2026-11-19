# Two additions to the label-sensitivity taus (reads saved scores only; nothing is rescored).
# 1. Weighted tau (scipy.stats.weightedtau, hyperbolic weights, rank = None: both rankings drive the weights)
#    between the EBCC-label and best-labeller method rankings, bio and overall, every tissue x dataset:
#    a swap among the top methods costs more than one lower down.
# 2. Sampling-noise baseline: the analysis/bootstrap_scores.py repeats (dissected dataset, 20,000-cell subsamples
#    stratified by species x EBCC label, EBCC labels). Kendall's tau between the method rankings of every pair
#    of repeats = how much rankings move from cell sampling alone. Compared with the label-change tau on the
#    same dataset, on the same methods, with the full six-metric bio score (bootstrap_bio6/summary/1_long.csv).
# Usage: python analysis/tau_extra.py  -> figures/cross_tissue/tau_weighted.csv, tau_noise_baseline.csv, tau_noise_summary.csv
import itertools, glob
import numpy as np, pandas as pd
from scipy.stats import kendalltau, weightedtau

import os
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = f'{ROOT}/figures/cross_tissue'
BATCH = ['PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn']
BIO = ['NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'isolated_label_F1', 'isolated_label_silhouette', 'cLISI']
BIO_BOOT = BIO   # all six: the raw bootstrap scores hold the isolated-label metrics (analysis/bootstrap_bio6_summarise.py)
TISSUES = {'bone_marrow': 'best_annotator', 'fetal_liver': 'best_annotator', 'breast_cancer': 'best_annotator'}
DATASETS = ['full', 'dissected', 'genes2', 'genes4', 'genes6', 'genes8']

def score_dir(t, d):
    dn = d + '_mouse' if d in ('full', 'dissected') else d
    return f'{ROOT}/{t}/data/{dn}/' if t == 'breast_cancer' or d.startswith('genes') else f'{ROOT}/{t}/data/05_evaluated_score/{dn}/'

def read_metrics(folder):
    rows = []
    for f in sorted(glob.glob(f'{folder}05.*.csv')):
        df = pd.read_csv(f)
        rows.append(pd.DataFrame({'metric': df.iloc[:, 0], 'value': pd.to_numeric(df.iloc[:, 1], errors='coerce'),
                                  'method': df.columns[1]}))
    return pd.concat(rows).dropna(subset=['value']) if rows else None

def scores(m, bio = BIO):
    g = lambda ms: m[m.metric.isin(ms)].groupby('method').value.mean()
    s = pd.DataFrame({'bio': g(bio), 'batch': g(BATCH)})
    s['overall'] = 0.6 * s.bio + 0.4 * s.batch
    return s

def wtau(x, y): return weightedtau(x, y).statistic

# one-sided permutation p for weighted tau (H1: tau_w > 0): shuffle y, count statistics >= observed
def wtau_perm_p(x, y, n = 9999, seed = 1):
    rng = np.random.default_rng(seed)
    obs = wtau(x, y); y = np.asarray(y)
    hits = sum(wtau(x, rng.permutation(y)) >= obs for _ in range(n))
    return (hits + 1) / (n + 1)

# ---- 1. weighted tau, label change, every tissue x dataset ----
out = []
for t, ls in TISSUES.items():
    for d in DATASETS:
        e, a = read_metrics(score_dir(t, d)), read_metrics(f'{score_dir(t, d)}label_{ls}/')
        if e is None or a is None: continue
        j = scores(e).join(scores(a), lsuffix = '_ebcc', rsuffix = '_alt', how = 'inner').dropna()
        for sc in ('bio', 'overall'):
            out.append(dict(tissue = t, dataset = d, score = sc, n_methods = len(j),
                            tau = kendalltau(j[f'{sc}_ebcc'], j[f'{sc}_alt']).statistic,
                            weighted_tau = wtau(j[f'{sc}_ebcc'], j[f'{sc}_alt']),
                            weighted_tau_p = wtau_perm_p(j[f'{sc}_ebcc'].to_numpy(), j[f'{sc}_alt'].to_numpy())))
w = pd.DataFrame(out)
w.to_csv(f'{FIG}/tau_weighted.csv', index = False)
print(w.round(3).to_string(index = False))

# ---- 2. sampling-noise baseline from the bootstrap repeats (dissected dataset) ----
boot = pd.read_csv(f'{ROOT}/bootstrap_bio6/summary/1_long.csv')
boot = boot[boot.repeat.astype(str) != 'full']
noise, label_same_basis = [], []
for t in ['colon', 'bone_marrow', 'fetal_liver', 'breast_cancer']:
    b = boot[boot.dataset == t]
    reps = {r: scores(g.rename(columns = {'score': 'value'}), BIO_BOOT) for r, g in b.groupby('repeat')}
    keys = sorted(reps)
    for r1, r2 in itertools.combinations(keys, 2):
        j = reps[r1].join(reps[r2], lsuffix = '_1', rsuffix = '_2', how = 'inner').dropna()
        for sc in ('bio', 'overall'):
            noise.append(dict(tissue = t, rep1 = r1, rep2 = r2, score = sc, n_methods = len(j),
                              tau = kendalltau(j[f'{sc}_1'], j[f'{sc}_2']).statistic,
                              weighted_tau = wtau(j[f'{sc}_1'], j[f'{sc}_2'])))
    if t in TISSUES:   # label change on the dissected dataset, bootstrap metric set and methods
        methods = set(b.method)
        e = read_metrics(score_dir(t, 'dissected')); a = read_metrics(f"{score_dir(t, 'dissected')}label_{TISSUES[t]}/")
        j = scores(e, BIO_BOOT).join(scores(a, BIO_BOOT), lsuffix = '_ebcc', rsuffix = '_alt', how = 'inner').dropna()
        j = j[j.index.isin(methods)]
        for sc in ('bio', 'overall'):
            label_same_basis.append(dict(tissue = t, score = sc, n_methods_label = len(j),
                                         label_tau = kendalltau(j[f'{sc}_ebcc'], j[f'{sc}_alt']).statistic,
                                         label_weighted_tau = wtau(j[f'{sc}_ebcc'], j[f'{sc}_alt'])))
noise = pd.DataFrame(noise)
noise.to_csv(f'{FIG}/tau_noise_baseline.csv', index = False)
q = lambda p: (lambda x: x.quantile(p))
summ = (noise.groupby(['tissue', 'score'])
        .agg(n_pairs = ('tau', 'size'), n_methods = ('n_methods', 'median'),
             noise_tau_median = ('tau', 'median'), noise_tau_q05 = ('tau', q(0.05)), noise_tau_q95 = ('tau', q(0.95)),
             noise_wtau_median = ('weighted_tau', 'median'), noise_wtau_q05 = ('weighted_tau', q(0.05)),
             noise_wtau_q95 = ('weighted_tau', q(0.95)))
        .reset_index()
        .merge(pd.DataFrame(label_same_basis), on = ['tissue', 'score'], how = 'left'))
# where the label-change tau sits in the noise distribution (share of repeat pairs with a lower tau)
summ['label_tau_percentile'] = [np.nan if np.isnan(r.label_tau) else
                                (noise[(noise.tissue == r.tissue) & (noise.score == r.score)].tau < r.label_tau).mean()
                                for r in summ.itertuples()]
summ.to_csv(f'{FIG}/tau_noise_summary.csv', index = False)
print(summ.round(3).to_string(index = False))
