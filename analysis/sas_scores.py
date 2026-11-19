# Seurat alignment score (SAS) of every integration, all datasets except full tissue.
#   python analysis/sas_scores.py run          # all jobs, JOBS at a time, then the summary table
#   python analysis/sas_scores.py job TISSUE DATASET SCRIPT
#   python analysis/sas_scores.py summarise
# SAS (Butler et al. 2018), as in the scGLUE Python implementation (scglue.metrics.seurat_alignment_score):
# the species are subsampled to the size of the smaller one, K = 1% of the subsampled cells,
# x = mean number of same-species cells among each cell's K nearest neighbours (self excluded), and
#   SAS = 1 - (x - K/N) / (K - K/N),  N = number of species; capped at 1, averaged over N_REPEATS draws.
# Each job runs the method's own 05 evaluation script with utils.score swapped for SAS, so the
# embedding, the evaluation cells and the input files are exactly those of the saved scib scores.
# Graph-only methods (output_type 'knn': BBKNN variants, SAMap) have no embedding: SAS is NA.
# Output: sas/<tissue>/<dataset>/<script>.csv, sas/sas_all.csv, log sas/run.log
import os, sys, time, glob, runpy, subprocess
import numpy as np
import pandas as pd

# ---------------- config ----------------
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/sas'
PY = os.path.join(os.environ.get('CONDA_ENVS', os.path.expanduser('~/miniconda3/envs')), 'scib/bin/python')
NEIGHBOR_FRAC = 0.01       # K = 1% of the subsampled cells
N_REPEATS = 4              # subsampling draws averaged (scGLUE default)
SEED = 20261030
SPECIES = 'orig.ident'     # species column (hsapiens / mmusculus)
JOBS = 6                   # jobs at once (the bootstrap runs alongside)
THREADS = 4                # threads per job
MIN_FREE_GB = 40           # wait to start a job until this much memory is available
DATASETS = {               # dataset names as each tissue's 05 scripts take them; full tissue left out
    'bone_marrow':   ['dissected', 'genes2', 'genes4', 'genes6', 'genes8'],
    'colon':         ['dissected_mouse', 'genes2', 'genes4', 'genes6', 'genes8'],
    'breast_cancer': ['dissected_mouse', 'genes2', 'genes4', 'genes6', 'genes8'],
    'fetal_liver':   ['dissected_mouse', 'genes2', 'genes4', 'genes6', 'genes8'],
}

def scripts(tissue):
    """every 05 evaluation script: the three unintegrated inputs and the 23 methods"""
    return sorted(os.path.basename(f) for f in glob.glob(f'{ROOT}/{tissue}/05.*.py') if 'scIB' not in f)

def out_file(tissue, dataset, script): return f'{OUT}/{tissue}/{dataset}/{script[:-3]}.csv'

# ---------------- the score ----------------
def seurat_alignment_score(x, y, neighbor_frac = NEIGHBOR_FRAC, n_repeats = N_REPEATS, seed = SEED, n_jobs = THREADS):
    from sklearn.neighbors import NearestNeighbors
    rs = np.random.RandomState(seed)
    idx_list = [np.where(y == u)[0] for u in np.unique(y)]
    n = len(idx_list)
    min_size = min(idx.size for idx in idx_list)
    out = []
    for _ in range(n_repeats):
        sub = np.concatenate([rs.choice(idx, min_size, replace = False) for idx in idx_list])
        xs, ys = x[sub], y[sub]
        k = max(round(sub.size * neighbor_frac), 1)
        nni = NearestNeighbors(n_neighbors = k + 1, n_jobs = n_jobs).fit(xs).kneighbors(xs, return_distance = False)
        same = (ys[nni[:, 1:]] == ys[:, None]).sum(axis = 1).mean()
        out.append(min((k - same) * n / (k * (n - 1)), 1))
    return float(np.mean(out)), k, min_size, n

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
            x = np.asarray(adata_int.obsm[embed], dtype = np.float32)
            y = adata_int.obs[SPECIES].astype(str).to_numpy()
            sas, k, min_size, n = seurat_alignment_score(x, y)
            res.update(sas = sas, k = k, min_species_cells = min_size, n_species = n)
        else:
            res.update(sas = np.nan)
        res['sas_seconds'] = time.time() - t0
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
    print(f'saved {out}: SAS {res["sas"]} ({res["seconds"]:.0f} s)')

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
    # smallest tissue first, so the first timings come in quickly
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

def summarise():
    files = glob.glob(f'{OUT}/*/*/*.csv')
    if not files:
        return
    d = pd.concat([pd.read_csv(f) for f in files]).sort_values(['tissue', 'dataset', 'script'])
    d.to_csv(f'{OUT}/sas_all.csv', index = False)
    log(f'summary: {len(d)} scores -> {OUT}/sas_all.csv')

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'job':
        job(*sys.argv[2:5])
    elif mode == 'run':
        run()
    elif mode == 'summarise':
        summarise()
