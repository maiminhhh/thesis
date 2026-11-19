# Extra integration metrics, alongside the scib scores (which are not changed):
#   MAP  mean average precision of cell-type labels among each cell's neighbours (bio)
#   NC   neighbour conservation: Jaccard overlap of each cell's neighbours in its own species'
#        PCA (before) and in the integrated embedding (after), within each species (bio, label-free)
#   LISI-CoM      harmonic mean of the saved scib iLISI and cLISI (combined)
#   Silhouette F1 harmonic mean of the saved scib ASW_label/batch and ASW_label (combined)
# MAP and NC follow the scGLUE implementation (scglue.metrics.mean_average_precision /
# neighbor_conservation, v0.3.2), the one STACAME uses: exact neighbours (sklearn), K = 1% of all
# cells (NC: 1% of all cells, searched within each species), MAP drops each cell's first neighbour.
# Neighbours are found in chunks so memory stays small; the values are those of the scGLUE code.
# "Before" for NC = each species' own PCA: utils.joint_pca's steps (data layer, scale, PCA) on one species.
# Each job runs the method's own 05 script with utils.score swapped, so embedding and cells are
# exactly those of the saved scib scores. Graph-only methods (no embedding): MAP / NC are NA.
#   python analysis/extra_metrics.py run          # MAP / NC jobs, then the summary
#   python analysis/extra_metrics.py job TISSUE DATASET SCRIPT
#   python analysis/extra_metrics.py summarise    # combined scores from the saved scib scores + the table
# Output: extra_metrics/<tissue>/<dataset>/<script>.csv, extra_metrics/extra_metrics_all.csv
import os, sys, time, glob, runpy, subprocess
import numpy as np
import pandas as pd

# ---------------- config ----------------
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/extra_metrics'
PY = os.path.join(os.environ.get('CONDA_ENVS', os.path.expanduser('~/miniconda3/envs')), 'scib/bin/python')
NEIGHBOR_FRAC = 0.01       # K = 1% of the cells (scGLUE default)
CHUNK = 4000               # query cells per neighbour search
SPECIES = 'orig.ident'
LABEL = 'ebcc_label'
JOBS = 10                  # jobs at once
THREADS = 4
MIN_FREE_GB = 60
GENE_LEVELS = [f'genes{g}' for g in (2, 4, 6, 8)]
DATASETS = {               # MAP / NC: dissected and the reduced-gene datasets
    'bone_marrow':   ['dissected', *GENE_LEVELS],
    'colon':         ['dissected_mouse', *GENE_LEVELS],
    'breast_cancer': ['dissected_mouse', *GENE_LEVELS],
    'fetal_liver':   ['dissected_mouse', *GENE_LEVELS],
}
# saved scib scores, for the combined scores (every dataset but full tissue)
SCORE_DIRS = {
    'bone_marrow':   {'dissected': 'data/05_evaluated_score/dissected_mouse/', **{f'genes{g}': f'data/genes{g}/' for g in (2, 4, 6, 8)}},
    'colon':         {'dissected_mouse': 'data/05_evaluated_score/dissected_mouse/', **{f'genes{g}': f'data/genes{g}/' for g in (2, 4, 6, 8)}},
    'breast_cancer': {'dissected_mouse': 'data/dissected_mouse/', **{f'genes{g}': f'data/genes{g}/' for g in (2, 4, 6, 8)}},
    'fetal_liver':   {'dissected_mouse': 'data/05_evaluated_score/dissected_mouse/', **{f'genes{g}': f'data/genes{g}/' for g in (2, 4, 6, 8)}},
}

def scripts(tissue):
    return sorted(os.path.basename(f) for f in glob.glob(f'{ROOT}/{tissue}/05.*.py') if 'scIB' not in f)

def out_file(tissue, dataset, script): return f'{OUT}/{tissue}/{dataset}/{script[:-3]}.csv'

# ---------------- metrics ----------------
def _knn(x, k):
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors = min(x.shape[0], k + 1), n_jobs = THREADS).fit(x)
    for s in range(0, x.shape[0], CHUNK):
        yield s, nn.kneighbors(x[s:s + CHUNK], return_distance = False)

def mean_average_precision(x, y):
    k = max(round(y.shape[0] * NEIGHBOR_FRAC), 1)
    ap = np.empty(y.shape[0])
    for s, nni in _knn(x, k):
        match = y[nni[:, 1:]] == y[s:s + nni.shape[0], None]
        prec = np.cumsum(match, axis = 1) / np.arange(1, match.shape[1] + 1)
        n = match.sum(axis = 1)
        ap[s:s + nni.shape[0]] = np.where(n > 0, (prec * match).sum(axis = 1) / np.maximum(n, 1), 0.0)
    return float(ap.mean()), k

def neighbor_conservation(x, y, batch):
    k = max(round(x.shape[0] * NEIGHBOR_FRAC), 1)
    per_batch = []
    for b in np.unique(batch):
        m = batch == b
        x_, y_ = x[m], y[m]
        jac = np.empty(m.sum())
        for (s, a), (_, c) in zip(_knn(x_, k), _knn(y_, k)):
            rows = np.arange(s, s + a.shape[0])[:, None]
            a = np.where(a == rows, -1, a); c = np.where(c == rows, -2, c)   # remove self (as setdiag(0))
            na, nc = (a >= 0).sum(1), (c >= 0).sum(1)
            both = np.sort(np.concatenate([a, c], axis = 1), axis = 1)
            inter = ((both[:, 1:] == both[:, :-1]) & (both[:, 1:] >= 0)).sum(1)
            jac[s:s + a.shape[0]] = inter / (na + nc - inter)
        per_batch.append(jac.mean())
    return float(np.mean(per_batch)), k

def species_pca(raw, species):
    """each species' own PCA, with utils.joint_pca's steps"""
    import scanpy as sc
    out = {}
    for b in np.unique(species):
        a = raw[species == b].copy()
        layer = next((l for l in a.layers if l.startswith('data')), None)
        if layer is not None:
            a.X = a.layers[layer].copy()
        sc.pp.scale(a, max_value = 10)
        sc.pp.pca(a)
        out[b] = a.obsm['X_pca'].astype(np.float32)
    return out

class _Captured(Exception):
    pass

def job(tissue, dataset, script):
    os.chdir(f'{ROOT}/{tissue}')
    sys.path.insert(0, os.getcwd())
    sys.argv = [script, dataset]
    import utils
    captured = []
    def score(adata_int, method, output_type, raw, embed = None, **kw):
        t0 = time.time()
        res = dict(tissue = tissue, dataset = dataset, script = script, method = method, output_type = output_type,
                   embed = embed, n_cells = adata_int.n_obs)
        if output_type == 'embed':
            assert (raw.obs_names == adata_int.obs_names).all(), 'raw and integrated cells differ'
            x = np.asarray(adata_int.obsm[embed], dtype = np.float32)
            obs = adata_int.obs if LABEL in adata_int.obs else raw.obs   # some integrated objects carry no labels
            lab = obs[LABEL].astype(str).to_numpy()
            sp = adata_int.obs[SPECIES].astype(str).to_numpy()
            res['MAP'], res['k_MAP'] = mean_average_precision(x, lab)
            res['MAP_seconds'] = time.time() - t0
            # NC: before = own-species PCA (rows in the same order as x within each species)
            pcs = species_pca(raw, sp)
            before = np.zeros((x.shape[0], max(p.shape[1] for p in pcs.values())), dtype = np.float32)
            for b, p in pcs.items():
                before[sp == b, :p.shape[1]] = p
            res['NC'], res['k_NC'] = neighbor_conservation(before, x, sp)
        else:
            res.update(MAP = np.nan, NC = np.nan)
        res['metric_seconds'] = time.time() - t0
        captured.append(res)
        raise _Captured
    utils.score = score
    t0 = time.time()
    try:
        runpy.run_path(script, run_name = '__main__')
    except _Captured:
        pass
    assert captured, f'{script}: utils.score was never called'
    res = dict(captured[0], seconds = time.time() - t0)
    out = out_file(tissue, dataset, script)
    os.makedirs(os.path.dirname(out), exist_ok = True)
    pd.DataFrame([res]).to_csv(out + '.tmp', index = False); os.replace(out + '.tmp', out)
    print(f'saved {out}: MAP {res["MAP"]}, NC {res["NC"]} ({res["seconds"]:.0f} s)')

# ---------------- runner ----------------
def free_gb():
    for line in open('/proc/meminfo'):
        if line.startswith('MemAvailable'):
            return int(line.split()[1]) / 1e6
    return 0

def log(msg):
    t = time.strftime("%H:%M:%S")   # '=====' markers start the line, for ~/bin/watch_job.sh
    line = f'{msg} ({t})' if msg.startswith('=====') else f'{t} {msg}'
    print(line, flush = True)
    with open(f'{OUT}/run.log', 'a') as f:
        f.write(line + '\n')

def run():
    os.makedirs(OUT, exist_ok = True)
    jobs = [(t, d, s) for t in DATASETS for d in DATASETS[t] for s in scripts(t)]
    todo = [j for j in jobs if not os.path.exists(out_file(*j))]
    log(f'===== start: {len(todo)} jobs ({len(jobs) - len(todo)} already done)')
    env = dict(os.environ, **{v: str(THREADS) for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                                        'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS', 'N_THREADS', 'N_CPU')})
    running, failed = {}, []
    while todo or running:
        for p, j in list(running.items()):
            if p.poll() is not None:
                del running[p]
                if p.returncode != 0:
                    failed.append(j)
                log(f'{"done  " if p.returncode == 0 else "FAILED"} {" ".join(j)}')
        while todo and len(running) < JOBS and free_gb() > MIN_FREE_GB:
            j = todo.pop(0)
            lf = f'{OUT}/logs/{j[0]}/{j[1]}/{j[2][:-3]}.log'
            os.makedirs(os.path.dirname(lf), exist_ok = True)
            p = subprocess.Popen([PY, os.path.abspath(__file__), 'job', *j], env = env,
                                 stdout = open(lf, 'w'), stderr = subprocess.STDOUT, cwd = ROOT)
            running[p] = j
            log(f'start  {" ".join(j)}')
        time.sleep(10)
    summarise()
    if failed:
        log(f'===== FAILED: {len(failed)} jobs: {failed}')
        sys.exit(1)
    log('===== ALL DONE')

def f1(a, b):
    return 2 * a * b / (a + b)

def summarise():
    os.makedirs(OUT, exist_ok = True)
    rows = []
    for t, dirs in SCORE_DIRS.items():
        for d, sd in dirs.items():
            for f in glob.glob(f'{ROOT}/{t}/{sd}05.*.csv'):
                s = pd.read_csv(f, index_col = 0)
                v = s.iloc[:, 0]
                g = lambda m: float(v[m]) if m in v.index and pd.notna(v[m]) else np.nan
                rows.append(dict(tissue = t, dataset = d, method = s.columns[0],
                                 iLISI = g('iLISI'), cLISI = g('cLISI'),
                                 ASW_batch = g('ASW_label/batch'), ASW_label = g('ASW_label')))
    comb = pd.DataFrame(rows)
    comb['LISI_CoM'] = f1(comb.iLISI, comb.cLISI)
    comb['silhouette_F1'] = f1(comb.ASW_batch, comb.ASW_label)
    files = glob.glob(f'{OUT}/*/*/*.csv')
    if files:
        nb = pd.concat([pd.read_csv(f) for f in files])[['tissue', 'dataset', 'method', 'MAP', 'NC', 'k_MAP', 'k_NC']]
        comb = comb.merge(nb, on = ['tissue', 'dataset', 'method'], how = 'left')
    comb.sort_values(['tissue', 'dataset', 'method']).to_csv(f'{OUT}/extra_metrics_all.csv', index = False)
    log(f'summary: {len(comb)} rows -> {OUT}/extra_metrics_all.csv')

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'job':
        job(*sys.argv[2:5])
    elif mode == 'run':
        run()
    elif mode == 'summarise':
        summarise()
