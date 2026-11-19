# Uncertainty of the benchmark scores and rankings by paired subsampling of cells (no integration rerun).
# Each repeat draws one subsample per dataset (fixed size, without replacement, stratified by species x
# EBCC label) and rescores every method and the unintegrated baseline on those same cells, with each
# tissue's own 05 script and utils.score unchanged: the runner swaps utils.eval_cells for the subsample,
# captures the result of utils.score and stops the script before it writes anything. Graph-only methods
# (BBKNN on PCA / UINMF / scVI latent, SAMap) are scored on their stored graph restricted to the subsample
# (utils.align subsets obsp). The clustering resolution for NMI/ARI is re-optimised in every repeat by
# scib, as in the pipeline. Overall = 0.6 bio + 0.4 batch on the metrics below, no min-max scaling.
#
#   python analysis/bootstrap_scores.py plan                 draw and save the subsamples (seeded)
#   python analysis/bootstrap_scores.py run [B]              check -> repeat 1 (timing) -> repeats 2..B, resumable
#   python analysis/bootstrap_scores.py job TISSUE REP SCRIPT   one scoring job (REP = full or a repeat number)
#   python analysis/bootstrap_scores.py summarise            outputs 1-6 from whatever repeats are complete
# Output: bootstrap/ (subsamples, raw scores, logs, summaries)
import os, sys, time, json, glob, runpy, subprocess
import numpy as np
import pandas as pd

# ---------------- config ----------------
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/bootstrap'
PY = os.path.join(os.environ.get('CONDA_ENVS', os.path.expanduser('~/miniconda3/envs')), 'scib/bin/python')
N_CELLS = 20000            # cells per subsample (per dataset)
MIN_PER_STRATUM = 5        # rare species x cell-type strata keep at least this many cells (or all they have)
B_DEFAULT = 50             # repeats
SEED_BASE = 20261030       # seed of repeat b in tissue t: SEED_BASE + 1000 * t + b
JOBS = 24                  # scoring jobs at once
THREADS = 3                # threads per job
MIN_FREE_GB = 30           # wait to start a job until this much memory is available
CHECK_TOL = 1e-4           # full-data rescoring must match the saved scores within this
LISI_TOL = 0.01            # except graph LISI, which varies between identical runs (seen up to 0.009)
BATCH_METRICS = ['PCR_batch', 'ASW_label/batch', 'iLISI', 'graph_conn']
BIO_METRICS = ['NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'cLISI']
W_BIO, W_BATCH = 0.6, 0.4
BASELINE = 'Unintegrated'
TISSUES = {   # dissected dataset of each tissue
    'bone_marrow':   dict(dataset = 'dissected',
                          concat = 'data/03_concatenated_data/03_concatenated_data_subset.h5ad',
                          eval_cells = 'data/04_integrated_data/dissected_mouse/eval_cells.txt',
                          scores = 'data/05_evaluated_score/dissected_mouse/'),
    'colon':         dict(dataset = 'dissected_mouse',
                          concat = 'data/03_concatenated_data/03_concatenated_data_subset.h5ad',
                          eval_cells = 'data/04_integrated_data/dissected_mouse/eval_cells.txt',
                          scores = 'data/05_evaluated_score/dissected_mouse/'),
    'breast_cancer': dict(dataset = 'dissected_mouse',
                          concat = 'data/dissected_mouse/03_concatenated_data.h5ad',
                          eval_cells = 'data/dissected_mouse/eval_cells.txt',
                          scores = 'data/dissected_mouse/'),
    'fetal_liver':   dict(dataset = 'dissected_mouse',
                          concat = 'data/03_concatenated_data/03_concatenated_data_subset.h5ad',
                          eval_cells = 'data/04_integrated_data/dissected_mouse/eval_cells.txt',
                          scores = 'data/05_evaluated_score/dissected_mouse/')}
CHECK_TISSUE = 'bone_marrow'   # the only dataset small enough to rescore at 100% here
# ----------------------------------------

def scripts(tissue):
    """the unintegrated baseline and the 23 methods (not the UINMF / scVI-latent unintegrated inputs)"""
    s = sorted(os.path.basename(f) for f in glob.glob(f'{ROOT}/{tissue}/05.*.py'))
    return [f for f in s if f.startswith('05.00_evaluate_raw.py') or
            not (f.startswith('05.00') or f.startswith('05.18') or 'scIB' in f)]

def sub_file(tissue, b): return f'{OUT}/{tissue}/subsamples/rep_{b:03d}.txt'
def score_file(tissue, rep, script): return f'{OUT}/{tissue}/scores/{rep}/{script[:-3]}.csv'

# ---------------- subsamples ----------------
def plan(B):
    import anndata as ad
    rows = []
    for t_i, (tissue, cfg) in enumerate(TISSUES.items()):
        os.makedirs(f'{OUT}/{tissue}/subsamples', exist_ok = True)
        cells = pd.Index(pd.read_csv(f'{ROOT}/{tissue}/{cfg["eval_cells"]}', header = None)[0].astype(str))
        obs = ad.read_h5ad(f'{ROOT}/{tissue}/{cfg["concat"]}', backed = 'r').obs
        cells = cells[cells.isin(obs.index)]
        strata = (obs.loc[cells, 'orig.ident'].astype(str) + '|' + obs.loc[cells, 'ebcc_label'].astype(str)).to_numpy()
        sizes = pd.Series(strata).value_counts()
        alloc = np.maximum(np.floor(sizes / sizes.sum() * N_CELLS), np.minimum(sizes, MIN_PER_STRATUM)).astype(int)
        for b in range(1, B + 1):
            f = sub_file(tissue, b)
            seed = SEED_BASE + 1000 * t_i + b
            if not os.path.exists(f):
                rng = np.random.default_rng(seed)
                pick = np.sort(np.concatenate([rng.choice(np.flatnonzero(strata == s), k, replace = False)
                                               for s, k in alloc.items()]))
                pd.Series(cells[pick]).to_csv(f, index = False, header = False)
            rows.append(dict(tissue = tissue, rep = b, seed = seed, n_cells = int(alloc.sum()),
                             n_strata = len(alloc), n_eval_cells = len(cells), file = f))
    pd.DataFrame(rows).to_csv(f'{OUT}/subsamples_manifest.csv', index = False)
    print(f'{B} subsamples per dataset -> {OUT}/subsamples_manifest.csv')

# ---------------- one scoring job ----------------
class _Captured(Exception):
    pass

def job(tissue, rep, script):
    cfg = TISSUES[tissue]
    os.chdir(f'{ROOT}/{tissue}')
    sys.path.insert(0, f'{ROOT}/{tissue}')
    sys.argv = [script, cfg['dataset']]
    import utils
    sub = None if rep == 'full' else pd.Index(pd.read_csv(sub_file(tissue, int(rep)), header = None)[0].astype(str))
    eval_cells0, score0, captured = utils.eval_cells, utils.score, []
    def eval_cells(*a, **k):
        cells = eval_cells0(*a, **k)
        return cells if sub is None else cells[cells.isin(sub)]
    def score(*a, **k):
        captured.append(score0(*a, **k))
        raise _Captured
    utils.eval_cells, utils.score = eval_cells, score
    t0 = time.time()
    try:
        runpy.run_path(script, run_name = '__main__')
    except _Captured:
        pass
    assert captured, f'{script}: utils.score was never called'
    out = score_file(tissue, rep, script)
    os.makedirs(os.path.dirname(out), exist_ok = True)
    res = captured[0]
    res.loc['_seconds'] = time.time() - t0
    res.loc['_n_cells'] = len(sub) if sub is not None else np.nan
    res.to_csv(out + '.tmp'); os.replace(out + '.tmp', out)
    print(f'saved {out} ({time.time() - t0:.0f} s)')

# ---------------- runner ----------------
def free_gb():
    for line in open('/proc/meminfo'):
        if line.startswith('MemAvailable'):
            return int(line.split()[1]) / 1e6
    return 0

def log(msg):
    line = f'{time.strftime("%H:%M:%S")} {msg}'
    print(line, flush = True)
    with open(f'{OUT}/run.log', 'a') as f:
        f.write(line + '\n')

def running_elsewhere():
    """pids of job processes left running by an earlier runner, by job"""
    out = subprocess.run(['pgrep', '-af', 'analysis/bootstrap_scores.py job'], capture_output = True, text = True).stdout
    found = {}
    for line in out.splitlines():
        f = line.split()
        if f[-4] == 'job':
            found[(f[-3], f[-2], f[-1])] = int(f[0])
    return found

def n_jobs():
    """JOBS, or the number in OUT/JOBS (read each loop, so it can change while running)"""
    try:
        return int(open(f'{OUT}/JOBS').read())
    except (OSError, ValueError):
        return JOBS

def pool(jobs, label, every = None):
    """run (tissue, rep, script) jobs, n_jobs() at a time; skips finished ones and ones an earlier
    runner left running (they count against the slots); calls every() every 20 min; returns the failures"""
    orphans = {j: p for j, p in running_elsewhere().items()}
    todo = [j for j in jobs if not os.path.exists(score_file(*j)) and tuple(map(str, j)) not in orphans]
    log(f'===== {label}: {len(todo)} jobs ({len(jobs) - len(todo)} already done or running)')
    last = time.time()
    env = dict(os.environ, **{v: str(THREADS) for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                                        'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS', 'N_THREADS', 'N_CPU')})
    running, failed = {}, []
    while todo or running:
        for p, j in list(running.items()):
            if p.poll() is not None:
                del running[p]
                (failed.append(j) if p.returncode != 0 else None)
                log(f'{"done  " if p.returncode == 0 else "FAILED"} {j[0]} {j[1]} {j[2]}')
        orphans = {j: p for j, p in orphans.items() if os.path.exists(f'/proc/{p}')}
        if every and time.time() - last > 1200:
            every(); last = time.time()
        while todo and len(running) + len(orphans) < n_jobs() and free_gb() > MIN_FREE_GB:
            j = todo.pop(0)
            lf = f'{OUT}/logs/{j[0]}/{j[1]}/{j[2][:-3]}.log'
            os.makedirs(os.path.dirname(lf), exist_ok = True)
            p = subprocess.Popen([PY, os.path.abspath(__file__), 'job', *map(str, j)], env = env,
                                 stdout = open(lf, 'w'), stderr = subprocess.STDOUT, cwd = ROOT)
            running[p] = j
            log(f'start  {j[0]} {j[1]} {j[2]}')
        time.sleep(10)
    return failed

def check_full():
    """rescoring at 100% of cells must reproduce the saved scores"""
    rows = []
    for s in scripts(CHECK_TISSUE):
        new = pd.read_csv(score_file(CHECK_TISSUE, 'full', s), index_col = 0)
        method = new.columns[0]
        old = None
        for f in glob.glob(f'{ROOT}/{CHECK_TISSUE}/{TISSUES[CHECK_TISSUE]["scores"]}05.*.csv'):
            d = pd.read_csv(f, index_col = 0)
            if d.columns[0] == method:
                old = d
        for m in BATCH_METRICS + BIO_METRICS:
            rows.append(dict(method = method, metric = m, saved = old.loc[m].iloc[0], rescored = new.loc[m].iloc[0]))
    d = pd.DataFrame(rows)
    d['diff'] = (d.rescored - d.saved).abs()
    d.to_csv(f'{OUT}/check_full_{CHECK_TISSUE}.csv', index = False)
    tol = np.where(d.metric.isin(['iLISI', 'cLISI']), LISI_TOL, CHECK_TOL)
    bad = d[(d['diff'] > tol) | (d.saved.isna() != d.rescored.isna())]
    return d, bad

def run(B):
    os.makedirs(OUT, exist_ok = True)
    plan(B)
    failed = pool([(CHECK_TISSUE, 'full', s) for s in scripts(CHECK_TISSUE)], f'check: {CHECK_TISSUE} at 100% of cells')
    if failed:
        log(f'===== FAILED: check jobs failed: {failed}'); sys.exit(1)
    d, bad = check_full()
    if len(bad):
        log(f'===== FAILED: full-data rescoring differs from the saved scores ({len(bad)} values, max diff '
            f'{bad["diff"].max():.3g}); see {OUT}/check_full_{CHECK_TISSUE}.csv'); sys.exit(1)
    log(f'===== check passed: {len(d)} values, max diff {d["diff"].max():.2g}')
    # one queue over all repeats (no wait for a repeat's slowest job), in repeat order, slow tissues
    # first within a repeat; the summary is refreshed every 20 min and at the end
    slow_first = ['fetal_liver', 'colon', 'breast_cancer', 'bone_marrow']
    jobs = [(t, b, s) for b in range(1, B + 1) for t in sorted(TISSUES, key = slow_first.index) for s in scripts(t)]
    all_failed = pool(jobs, f'repeats 1-{B}', every = summarise)
    summarise()
    if all_failed:
        log(f'===== FAILED: {len(all_failed)} jobs failed, e.g. {all_failed[:3]}'); sys.exit(1)
    log('===== ALL DONE')

# ---------------- summaries ----------------
def kendall_w(R):
    """R: datasets x methods rank matrix"""
    m, n = R.shape
    S = ((R.sum(0) - m * (n + 1) / 2) ** 2).sum()
    return 12 * S / (m ** 2 * (n ** 3 - n))

MAX_MISSING = 3   # most failed methods a repeat may have and still be used

def failed_jobs():
    """(tissue, repeat, script) of jobs whose last attempt failed, from the run log"""
    last = {}
    if os.path.exists(f'{OUT}/run.log'):
        for line in open(f'{OUT}/run.log'):
            f = line.split()
            if len(f) == 5 and f[1] in ('done', 'FAILED'):
                last[(f[2], f[3], f[4])] = f[1]
    return {j for j, v in last.items() if v == 'FAILED'}

def summarise():
    from scipy.stats import kendalltau
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = []
    for f in glob.glob(f'{OUT}/*/scores/*/*.csv'):
        tissue, _, rep = f.split('/')[-4:-1]
        d = pd.read_csv(f, index_col = 0)
        method = d.columns[0]
        for m in BATCH_METRICS + BIO_METRICS:
            rows.append(dict(dataset = tissue, repeat = rep, method = method, metric = m, score = d.loc[m].iloc[0]))
    if not rows:
        return
    long = pd.DataFrame(rows)
    # full-data point estimates: the saved primary scores
    full = []
    for tissue, cfg in TISSUES.items():
        for f in glob.glob(f'{ROOT}/{tissue}/{cfg["scores"]}05.*.csv'):
            d = pd.read_csv(f, index_col = 0)
            for m in BATCH_METRICS + BIO_METRICS:
                if m in d.index:
                    full.append(dict(dataset = tissue, repeat = 'point', method = d.columns[0], metric = m, score = d.loc[m].iloc[0]))
    full = pd.DataFrame(full)
    allm = pd.concat([long[long.repeat != 'full'], full])

    def components(df):
        w = df.pivot_table(index = ['dataset', 'repeat', 'method'], columns = 'metric', values = 'score')
        out = pd.DataFrame({'batch': w[BATCH_METRICS].mean(axis = 1), 'bio': w[BIO_METRICS].mean(axis = 1)})
        out['overall'] = W_BIO * out.bio + W_BATCH * out.batch
        return out.reset_index().melt(id_vars = ['dataset', 'repeat', 'method'], var_name = 'metric', value_name = 'score')
    comp = components(allm)
    tab = pd.concat([allm, comp])
    # repeats where a job failed (graph LISI on some subsampled BBKNN graphs) keep the other methods:
    # the failed method is missing from that repeat only, and ranks are among the methods scored.
    # A repeat is used once the baseline and all but at most MAX_MISSING methods are scored.
    n_methods = {t: len(scripts(t)) for t in TISSUES}
    counts = tab[tab.metric == 'overall'].groupby(['dataset', 'repeat']).method.nunique()
    has_base = tab[tab.method == BASELINE].groupby(['dataset', 'repeat']).size().index
    usable = counts[[n >= n_methods[t] - MAX_MISSING for (t, r), n in counts.items()]].index.intersection(has_base)
    usable = usable[[r == 'point' or os.path.exists(f'{OUT}/{t}/scores/{r}') and
                     all(os.path.exists(score_file(t, r, s)) or (t, r, s) in failed_jobs() for s in scripts(t))
                     for t, r in usable]]
    tab = tab.set_index(['dataset', 'repeat']).loc[lambda x: x.index.isin(usable)].reset_index()
    base = tab[tab.method == BASELINE][['dataset', 'repeat', 'metric', 'score']].rename(columns = {'score': 'base'})
    tab = tab.merge(base, on = ['dataset', 'repeat', 'metric'], how = 'left')
    tab['delta_from_unintegrated'] = tab.score - tab.base
    meth = tab[tab.method != BASELINE].copy()
    meth['rank'] = meth.groupby(['dataset', 'repeat', 'metric']).score.rank(ascending = False)
    os.makedirs(f'{OUT}/summary', exist_ok = True)
    meth.drop(columns = 'base').to_csv(f'{OUT}/summary/1_long.csv', index = False)
    reps = meth[meth.repeat != 'point']
    log_n = reps.groupby('dataset').repeat.nunique().to_dict()

    # 2. point estimate and 2.5 / 97.5% percentiles
    s = (reps[reps.metric.isin(['batch', 'bio', 'overall'])].groupby(['dataset', 'method', 'metric'])
         .agg(q025 = ('score', lambda x: x.quantile(0.025)), median = ('score', 'median'),
              q975 = ('score', lambda x: x.quantile(0.975)),
              delta_q025 = ('delta_from_unintegrated', lambda x: x.quantile(0.025)),
              delta_q975 = ('delta_from_unintegrated', lambda x: x.quantile(0.975)),
              rank_median = ('rank', 'median'), n_repeats = ('score', 'size')).reset_index())
    pt = meth[(meth.repeat == 'point') & meth.metric.isin(['batch', 'bio', 'overall'])][['dataset', 'method', 'metric', 'score', 'rank']]
    s = s.merge(pt.rename(columns = {'score': 'full_data', 'rank': 'full_rank'}), on = ['dataset', 'method', 'metric'], how = 'left')
    s.to_csv(f'{OUT}/summary/2_score_intervals.csv', index = False)
    pd.DataFrame(sorted(failed_jobs()), columns = ['dataset', 'repeat', 'script']).to_csv(
        f'{OUT}/summary/failed_jobs.csv', index = False)

    # 3. rank frequencies and 4. P(A beats B), overall score
    for t in sorted(reps.dataset.unique()):
        o = reps[(reps.dataset == t) & (reps.metric == 'overall')]
        order = o.groupby('method')['rank'].mean().sort_values().index
        rf = (pd.crosstab(o.method, o['rank'].round().astype(int), normalize = 'index').reindex(order))
        rf.to_csv(f'{OUT}/summary/3_rank_frequency_{t}.csv')
        fig, ax = plt.subplots(figsize = (11, 8))
        im = ax.imshow(rf.values, cmap = 'Blues', vmin = 0, vmax = 1, aspect = 'auto')
        ax.set_xticks(range(rf.shape[1])); ax.set_xticklabels(rf.columns)
        ax.set_yticks(range(rf.shape[0])); ax.set_yticklabels(rf.index)
        ax.set_xlabel('Rank by overall score (1 = best)'); ax.set_title(f'{t}: share of {o.repeat.nunique()} repeats at each rank')
        fig.colorbar(im, ax = ax, label = 'Share of repeats'); fig.tight_layout()
        fig.savefig(f'{OUT}/summary/3_rank_frequency_{t}.png', dpi = 200); plt.close(fig)
        w = o.pivot(index = 'repeat', columns = 'method', values = 'score')[order]
        P = pd.DataFrame({a: [(w[a] > w[b]).mean() for b in order] for a in order}, index = order).T
        P.index.name = 'P(row scores higher than column)'
        P.to_csv(f'{OUT}/summary/4_pairwise_{t}.csv')

    # 5. Kendall tau with published rankings (bootstrap/published/<paper>.csv: our_method, rank)
    pubs = {os.path.basename(f)[:-4]: pd.read_csv(f).dropna(subset = ['our_method']).set_index('our_method')['rank']
            for f in glob.glob(f'{OUT}/published/*.csv')}
    if pubs:
        o = reps[reps.metric == 'overall']
        res = []
        for name, pr in pubs.items():
            for t in sorted(o.dataset.unique()) + ['mean over datasets']:
                d = o if t == 'mean over datasets' else o[o.dataset == t]
                taus = []
                for r, g in d.groupby('repeat'):
                    ours = g.groupby('method')['rank'].mean()
                    shared = ours.index.intersection(pr.index)
                    taus.append(kendalltau(ours[shared], pr[shared])[0])
                res.append(dict(published = name, dataset = t, n_shared = len(shared), tau_median = np.median(taus),
                                tau_q025 = np.quantile(taus, 0.025), tau_q975 = np.quantile(taus, 0.975), n_repeats = len(taus)))
        names = list(pubs)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = pubs[names[i]], pubs[names[j]]
                shared = a.index.intersection(b.index)
                res.append(dict(published = f'{names[i]} vs {names[j]}', dataset = 'published', n_shared = len(shared),
                                tau_median = kendalltau(a[shared], b[shared])[0] if len(shared) > 2 else np.nan))
        pd.DataFrame(res).to_csv(f'{OUT}/summary/5_kendall_tau_published.csv', index = False)

    # 6. Kendall's W across datasets, per repeat number (repeats complete in every dataset)
    o = reps[reps.metric == 'overall']
    ws = []
    for r, g in o.groupby('repeat'):
        R = g.pivot(index = 'dataset', columns = 'method', values = 'score')
        if R.shape[0] == len(TISSUES):
            R = R.dropna(axis = 1).rank(axis = 1, ascending = False)   # methods scored in every dataset
            ws.append(dict(repeat = r, W = kendall_w(R.values), n_methods = R.shape[1]))
    if ws:
        w = pd.DataFrame(ws)
        pt_R = meth[(meth.repeat == 'point') & (meth.metric == 'overall')].pivot(index = 'dataset', columns = 'method', values = 'rank')
        w.to_csv(f'{OUT}/summary/6_kendall_w.csv', index = False)
        pd.DataFrame([dict(W_full_data = kendall_w(pt_R.dropna(axis = 1).values), W_median = w.W.median(),
                           W_q025 = w.W.quantile(0.025), W_q975 = w.W.quantile(0.975), n_repeats = len(w))]
                     ).to_csv(f'{OUT}/summary/6_kendall_w_summary.csv', index = False)
    print(f'summaries from {log_n} complete repeats -> {OUT}/summary/')

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'plan':
        plan(int(sys.argv[2]) if len(sys.argv) > 2 else B_DEFAULT)
    elif mode == 'job':
        job(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else B_DEFAULT)
    elif mode == 'summarise':
        summarise()
