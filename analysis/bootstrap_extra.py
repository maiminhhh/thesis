# Bootstrap of MAP, neighbour conservation (NC) and the Seurat alignment score (SAS), on the same paired
# subsamples as analysis/bootstrap_scores.py (bootstrap/<tissue>/subsamples/rep_XXX.txt), dissected tissue.
# Each job runs the method's 05 script with utils.eval_cells restricted to the subsample (as analysis/bootstrap_scores.py)
# and utils.score swapped for the metrics, computed exactly as in analysis/extra_metrics.py (MAP, NC) and analysis/sas_scores.py
# (SAS) on the subsampled cells: K = 1% of the subsample, NC "before" = each species' own PCA of the subsample.
# Embedding methods and the three unintegrated inputs only (graph-only methods have no embedding).
#   python analysis/bootstrap_extra.py run [B]                 all jobs, resumable, then the summary
#   python analysis/bootstrap_extra.py job TISSUE REP SCRIPT
#   python analysis/bootstrap_extra.py summarise
# Output: bootstrap_extra/<tissue>/<rep>/<script>.csv, bootstrap_extra/extra_long.csv, log bootstrap_extra/run.log
import os, sys, time, glob, runpy, subprocess
import numpy as np
import pandas as pd

# ---------------- config ----------------
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/bootstrap_extra'
PY = os.path.join(os.environ.get('CONDA_ENVS', os.path.expanduser('~/miniconda3/envs')), 'scib/bin/python')
B_DEFAULT = 50
JOBS = 64                  # default; the number in bootstrap_extra/JOBS overrides it, read each loop
THREADS = 2
MIN_FREE_GB = 40
DATASET = {'bone_marrow': 'dissected', 'colon': 'dissected_mouse', 'breast_cancer': 'dissected_mouse',
           'fetal_liver': 'dissected_mouse'}
GRAPH_ONLY = ('BBKNN', 'SAMap')

sys.path.insert(0, f'{ROOT}/analysis')
import extra_metrics as em, sas_scores as ss
em.THREADS = THREADS

def scripts(tissue):
    return sorted(os.path.basename(f) for f in glob.glob(f'{ROOT}/{tissue}/05.*.py')
                  if 'scIB' not in f and not any(g in os.path.basename(f) for g in GRAPH_ONLY))

def sub_file(tissue, b): return f'{ROOT}/bootstrap/{tissue}/subsamples/rep_{b:03d}.txt'
def out_file(tissue, rep, script): return f'{OUT}/{tissue}/{rep}/{script[:-3]}.csv'

class _Captured(Exception):
    pass

def job(tissue, rep, script):
    os.chdir(f'{ROOT}/{tissue}')
    sys.path.insert(0, f'{ROOT}/{tissue}')
    sys.argv = [script, DATASET[tissue]]
    import utils
    sub = pd.Index(pd.read_csv(sub_file(tissue, int(rep)), header = None)[0].astype(str))
    eval_cells0, captured = utils.eval_cells, []
    def eval_cells(*a, **k):
        cells = eval_cells0(*a, **k)
        return cells[cells.isin(sub)]
    def score(adata_int, method, output_type, raw, embed = None, **kw):
        t0 = time.time()
        assert output_type == 'embed', f'{script}: not an embedding method'
        assert (raw.obs_names == adata_int.obs_names).all(), 'raw and integrated cells differ'
        x = np.asarray(adata_int.obsm[embed], dtype = np.float32)
        obs = adata_int.obs if em.LABEL in adata_int.obs else raw.obs
        lab = obs[em.LABEL].astype(str).to_numpy()
        sp = adata_int.obs[em.SPECIES].astype(str).to_numpy()
        res = dict(tissue = tissue, rep = int(rep), script = script, method = method, n_cells = adata_int.n_obs)
        res['MAP'], res['k_MAP'] = em.mean_average_precision(x, lab)
        pcs = em.species_pca(raw, sp)
        before = np.zeros((x.shape[0], max(p.shape[1] for p in pcs.values())), dtype = np.float32)
        for b, p in pcs.items():
            before[sp == b, :p.shape[1]] = p
        res['NC'], res['k_NC'] = em.neighbor_conservation(before, x, sp)
        res['SAS'], res['k_SAS'], _, _ = ss.seurat_alignment_score(x, sp, n_jobs = THREADS)
        res['metric_seconds'] = time.time() - t0
        captured.append(res)
        raise _Captured
    utils.eval_cells, utils.score = eval_cells, score
    t0 = time.time()
    try:
        runpy.run_path(script, run_name = '__main__')
    except _Captured:
        pass
    assert captured, f'{script}: utils.score was never called'
    res = dict(captured[0], seconds = time.time() - t0)
    out = out_file(tissue, rep, script)
    os.makedirs(os.path.dirname(out), exist_ok = True)
    pd.DataFrame([res]).to_csv(out + '.tmp', index = False); os.replace(out + '.tmp', out)
    print(f'saved {out}: MAP {res["MAP"]:.3f} NC {res["NC"]:.3f} SAS {res["SAS"]:.3f} ({res["seconds"]:.0f} s)')

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
    slow_first = ['fetal_liver', 'colon', 'breast_cancer', 'bone_marrow']
    jobs = [(t, str(b), s) for b in range(1, B + 1) for t in slow_first for s in scripts(t)]
    todo = [j for j in jobs if not os.path.exists(out_file(*j))]
    log(f'===== start: {len(todo)} jobs ({len(jobs) - len(todo)} already done)')
    env = dict(os.environ, **{v: str(THREADS) for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                                        'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS', 'N_THREADS', 'N_CPU')})
    running, failed, last = {}, [], time.time()
    while todo or running:
        for p, j in list(running.items()):
            if p.poll() is not None:
                del running[p]
                if p.returncode != 0:
                    failed.append(j)
                log(f'{"done  " if p.returncode == 0 else "FAILED"} {" ".join(j)}')
        while todo and len(running) < n_jobs() and free_gb() > MIN_FREE_GB:
            j = todo.pop(0)
            lf = f'{OUT}/logs/{j[0]}/{j[1]}/{j[2][:-3]}.log'
            os.makedirs(os.path.dirname(lf), exist_ok = True)
            p = subprocess.Popen([PY, os.path.abspath(__file__), 'job', *j], env = env,
                                 stdout = open(lf, 'w'), stderr = subprocess.STDOUT, cwd = ROOT)
            running[p] = j
            log(f'start  {" ".join(j)}')
        if time.time() - last > 1200:
            summarise(); last = time.time()
        time.sleep(5)
    summarise()
    if failed:
        log(f'===== FAILED: {len(failed)} jobs, e.g. {failed[:3]}')
        sys.exit(1)
    log('===== ALL DONE')

def summarise():
    files = glob.glob(f'{OUT}/*/*/*.csv')
    if not files:
        return
    d = pd.concat([pd.read_csv(f) for f in files]).sort_values(['tissue', 'rep', 'method'])
    d.to_csv(f'{OUT}/extra_long.csv', index = False)
    log(f'summary: {len(d)} rows -> {OUT}/extra_long.csv')

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'job':
        job(*sys.argv[2:5])
    elif mode == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else B_DEFAULT)
    elif mode == 'summarise':
        summarise()
