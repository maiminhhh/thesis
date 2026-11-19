# Bootstrap of the reduced-gene datasets (genes2/4/6/8 = 20/40/60/80% of shared orthologues removed), as the
# dissected bootstrap (analysis/bootstrap_scores.py, analysis/bootstrap_extra.py): paired stratified subsamples of 20,000 cells
# (species x EBCC label), every method rescored on the same cells with its own 05 script and utils.score unchanged.
# The reduced-gene datasets keep fewer cells than dissected (cells without counts on the remaining genes are
# dropped), so each dataset gets its own subsamples, drawn as in bootstrap_scores.plan.
# One job = one method on one subsample: the scib metrics (utils.score) and MAP, NC, SAS (as analysis/extra_metrics.py,
# analysis/sas_scores.py, analysis/bootstrap_extra.py) on the same embedding and cells. Graph-only methods (BBKNN, SAMap) are left
# out: their stored graph cannot be subsampled.
#   python analysis/bootstrap_genes.py plan [B]                     draw and save the subsamples (seeded)
#   python analysis/bootstrap_genes.py run [B]                      all jobs, resumable, then the summary
#   python analysis/bootstrap_genes.py job TISSUE DATASET REP SCRIPT
#   python analysis/bootstrap_genes.py batch T,D,R,S T,D,R,S ...   several jobs in one process (numba compiles once)
#   python analysis/bootstrap_genes.py summarise
# Output: bootstrap_genes/<tissue>/<dataset>/{subsamples,scores}/, bootstrap_genes/genes_long.csv, run.log
import os, sys, time, glob, runpy, subprocess
import numpy as np
import pandas as pd

# ---------------- config ----------------
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/bootstrap_genes'
PY = os.path.join(os.environ.get('CONDA_ENVS', os.path.expanduser('~/miniconda3/envs')), 'scib/bin/python')
LEVELS = ['genes2', 'genes4', 'genes6', 'genes8']
TISSUES = ['bone_marrow', 'colon', 'breast_cancer', 'fetal_liver']
N_CELLS = 20000            # as analysis/bootstrap_scores.py
MIN_PER_STRATUM = 5
B_DEFAULT = 50
SEED_BASE = 20261101       # seed = SEED_BASE + 10000 * tissue index + 1000 * level index + repeat
JOBS = 32                  # default; the number in bootstrap_genes/JOBS overrides it, read each loop
THREADS = 2
MIN_FREE_GB = 40
GRAPH_ONLY = ('BBKNN', 'SAMap')
BATCH = 8                  # jobs per worker process: scanpy's neighbour search (pynndescent) spends ~150 s compiling
                           # with numba in every new process, so a worker runs several jobs; the scores are unchanged
# -----------------------------------------

sys.path.insert(0, f'{ROOT}/analysis')
import extra_metrics as em, sas_scores as ss
em.THREADS = THREADS

def paths(tissue, dataset):
    d = f'{ROOT}/{tissue}/data/{dataset}'
    return dict(concat = f'{d}/03_concatenated_data.h5ad', eval_cells = f'{d}/eval_cells.txt')

def scripts(tissue):
    """the unintegrated baseline and the embedding methods (as bootstrap_scores.scripts, without graph-only)"""
    s = sorted(os.path.basename(f) for f in glob.glob(f'{ROOT}/{tissue}/05.*.py'))
    s = [f for f in s if f.startswith('05.00_evaluate_raw.py') or
         not (f.startswith('05.00') or f.startswith('05.18') or 'scIB' in f)]
    return [f for f in s if not any(g in f for g in GRAPH_ONLY)]

def sub_file(tissue, dataset, b): return f'{OUT}/{tissue}/{dataset}/subsamples/rep_{b:03d}.txt'
def score_file(tissue, dataset, rep, script): return f'{OUT}/{tissue}/{dataset}/scores/{rep}/{script[:-3]}.csv'

# ---------------- subsamples ----------------
def plan(B):
    import anndata as ad
    rows = []
    for t_i, tissue in enumerate(TISSUES):
        for l_i, dataset in enumerate(LEVELS):
            p = paths(tissue, dataset)
            os.makedirs(os.path.dirname(sub_file(tissue, dataset, 1)), exist_ok = True)
            cells = pd.Index(pd.read_csv(p['eval_cells'], header = None)[0].astype(str))
            obs = ad.read_h5ad(p['concat'], backed = 'r').obs
            cells = cells[cells.isin(obs.index)]
            strata = (obs.loc[cells, 'orig.ident'].astype(str) + '|' + obs.loc[cells, 'ebcc_label'].astype(str)).to_numpy()
            sizes = pd.Series(strata).value_counts()
            alloc = np.maximum(np.floor(sizes / sizes.sum() * N_CELLS), np.minimum(sizes, MIN_PER_STRATUM)).astype(int)
            for b in range(1, B + 1):
                f = sub_file(tissue, dataset, b)
                seed = SEED_BASE + 10000 * t_i + 1000 * l_i + b
                if not os.path.exists(f):
                    rng = np.random.default_rng(seed)
                    pick = np.sort(np.concatenate([rng.choice(np.flatnonzero(strata == s), k, replace = False)
                                                   for s, k in alloc.items()]))
                    pd.Series(cells[pick]).to_csv(f, index = False, header = False)
                rows.append(dict(tissue = tissue, dataset = dataset, rep = b, seed = seed, n_cells = int(alloc.sum()),
                                 n_strata = len(alloc), n_eval_cells = len(cells), file = f))
            print(f'{tissue} {dataset}: {len(cells)} eval cells, {int(alloc.sum())} per subsample', flush = True)
    pd.DataFrame(rows).to_csv(f'{OUT}/subsamples_manifest.csv', index = False)

# ---------------- one scoring job ----------------
class _Captured(Exception):
    pass

def job(tissue, dataset, rep, script):
    t_start = time.time()
    os.chdir(f'{ROOT}/{tissue}')
    sys.path.insert(0, f'{ROOT}/{tissue}')
    sys.argv = [script, dataset]
    import utils
    sub = pd.Index(pd.read_csv(sub_file(tissue, dataset, int(rep)), header = None)[0].astype(str))
    eval_cells0, score0, captured = utils.eval_cells, utils.score, []
    def eval_cells(*a, **k):
        cells = eval_cells0(*a, **k)
        return cells[cells.isin(sub)]
    def score(adata_int, method, output_type, raw, embed = None, **kw):
        t0 = time.time()
        assert output_type == 'embed', f'{script}: not an embedding method'
        assert (raw.obs_names == adata_int.obs_names).all(), 'raw and integrated cells differ'
        # MAP, NC, SAS first (as analysis/bootstrap_extra.py), on a copy of the embedding, before scib touches adata_int
        x = np.asarray(adata_int.obsm[embed], dtype = np.float32).copy()
        obs = adata_int.obs if em.LABEL in adata_int.obs else raw.obs
        lab = obs[em.LABEL].astype(str).to_numpy()
        sp = adata_int.obs[em.SPECIES].astype(str).to_numpy()
        extra = {}
        extra['MAP'], extra['k_MAP'] = em.mean_average_precision(x, lab)
        pcs = em.species_pca(raw, sp)
        before = np.zeros((x.shape[0], max(p.shape[1] for p in pcs.values())), dtype = np.float32)
        for b, p in pcs.items():
            before[sp == b, :p.shape[1]] = p
        extra['NC'], extra['k_NC'] = em.neighbor_conservation(before, x, sp)
        extra['SAS'], extra['k_SAS'], _, _ = ss.seurat_alignment_score(x, sp, n_jobs = THREADS)
        t1 = time.time()
        res = score0(adata_int = adata_int, method = method, output_type = output_type, raw = raw, embed = embed, **kw)
        for k, v in extra.items():
            res.loc[k] = v
        res.loc['_load_seconds'] = t0 - t_start
        res.loc['_extra_seconds'] = t1 - t0
        res.loc['_scib_seconds'] = time.time() - t1
        res.loc['_n_cells'] = adata_int.n_obs
        captured.append(res)
        raise _Captured
    utils.eval_cells, utils.score = eval_cells, score
    try:
        runpy.run_path(script, run_name = '__main__')
    except _Captured:
        pass
    assert captured, f'{script}: utils.score was never called'
    out = score_file(tissue, dataset, rep, script)
    os.makedirs(os.path.dirname(out), exist_ok = True)
    res = captured[0]
    res.loc['_seconds'] = time.time() - t_start
    res.to_csv(out + '.tmp'); os.replace(out + '.tmp', out)
    print(f'saved {out} ({res.loc["_seconds"].iloc[0]:.0f} s)')

def batch(specs):
    """run several jobs in this process; returns the failed ones"""
    import gc, traceback
    failed = []
    for spec in specs:
        j = spec.split(',')
        if os.path.exists(score_file(*j)):
            continue
        for m in ('utils', 'datasets'):          # each tissue has its own utils.py / datasets.py
            sys.modules.pop(m, None)
        sys.path[:] = [p for p in sys.path if not p.startswith(f'{ROOT}/') or p in (ROOT, f'{ROOT}/analysis')]
        try:
            job(*j)
            print(f'BATCH done {spec}', flush = True)
        except Exception:
            traceback.print_exc()
            print(f'BATCH FAILED {spec}', flush = True)
            failed.append(spec)
        os.chdir(ROOT)
        gc.collect()
    return failed

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

def n_jobs():
    try:
        return int(open(f'{OUT}/JOBS').read())
    except (OSError, ValueError):
        return JOBS

def run(B):
    os.makedirs(OUT, exist_ok = True)
    plan(B)
    slow_first = ['fetal_liver', 'colon', 'breast_cancer', 'bone_marrow']
    # repeat-major, so the first repeats of every dataset finish early
    jobs = [(t, d, str(b), s) for b in range(1, B + 1) for t in slow_first for d in LEVELS for s in scripts(t)]
    todo = [j for j in jobs if not os.path.exists(score_file(*j))]
    log(f'===== start: {len(todo)} jobs ({len(jobs) - len(todo)} already done)')
    env = dict(os.environ, **{v: str(THREADS) for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                                        'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS', 'N_THREADS', 'N_CPU')})
    running, failed, last = {}, [], time.time()
    while todo or running:
        for p, (lf, js) in list(running.items()):
            if p.poll() is not None:
                del running[p]
                bad = [j for j in js if not os.path.exists(score_file(*j))]
                failed += bad
                for j in js:
                    log(f'{"FAILED" if j in bad else "done  "} {" ".join(j)}')
        while todo and len(running) < n_jobs() and free_gb() > MIN_FREE_GB:
            js, todo = todo[:BATCH], todo[BATCH:]
            lf = f'{OUT}/logs/batch_{js[0][0]}_{js[0][1]}_{js[0][2]}_{js[0][3][:5]}_{int(time.time() * 1000) % 10**8}.log'
            os.makedirs(os.path.dirname(lf), exist_ok = True)
            p = subprocess.Popen([PY, os.path.abspath(__file__), 'batch', *[','.join(j) for j in js]], env = env,
                                 stdout = open(lf, 'w'), stderr = subprocess.STDOUT, cwd = ROOT)
            running[p] = (lf, js)
            log(f'start  batch of {len(js)}: {" ".join(js[0])} ...')
        if time.time() - last > 1800:
            summarise(); last = time.time()
        time.sleep(5)
    summarise()
    if failed:
        log(f'===== FAILED: {len(failed)} jobs, e.g. {failed[:3]}')
        sys.exit(1)
    log('===== ALL DONE')

def summarise():
    rows = []
    for f in glob.glob(f'{OUT}/*/*/scores/*/*.csv'):
        tissue, dataset, _, rep = f[len(OUT) + 1:].split('/')[:4]
        s = pd.read_csv(f, index_col = 0)
        v = s.iloc[:, 0]
        rows.append(pd.DataFrame(dict(tissue = tissue, dataset = dataset, rep = int(rep), method = s.columns[0],
                                      metric = v.index, score = pd.to_numeric(v.values, errors = 'coerce'))))
    if not rows:
        return
    d = pd.concat(rows).dropna(subset = ['score']).sort_values(['tissue', 'dataset', 'rep', 'method', 'metric'])
    d.to_csv(f'{OUT}/genes_long.csv', index = False)
    log(f'summary: {d.groupby(["tissue", "dataset", "rep", "method"]).ngroups} scored jobs -> {OUT}/genes_long.csv')

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'plan':
        os.makedirs(OUT, exist_ok = True); plan(int(sys.argv[2]) if len(sys.argv) > 2 else B_DEFAULT)
    elif mode == 'job':
        job(*sys.argv[2:6])
    elif mode == 'batch':
        sys.exit(1 if batch(sys.argv[2:]) else 0)
    elif mode == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else B_DEFAULT)
    elif mode == 'summarise':
        summarise()
