# MAP under the best-annotator labels, for the label-robustness check with the +MAP +SAS score (analysis/label_robustness.R,
# extra = TRUE). SAS is label-free (species mixing only), so its saved values (sas/sas_all.csv) serve both labellings.
# MAP exactly as analysis/extra_metrics.py (scGLUE mean_average_precision: exact neighbours, K = 1% of the cells, first neighbour
# dropped), on the same embedding and cells: each job runs the method's own 05 script with utils.score swapped. One
# neighbour search per job gives MAP under both labellings: EBCC (must equal the saved extra_metrics MAP: the check) and
# the best annotator (per species, the obs column of label_sensitivity_choice.json, built as utils.score does).
# Embedding methods and the three unintegrated inputs; graph-only methods (BBKNN, SAMap) have no MAP, as in the benchmark.
# Datasets with MAP / SAS in the benchmark: dissected and genes2-8 (no full tissue); colon has no alternative labelling.
#   python analysis/label_map_alt.py run           all jobs (resumable), then the summary
#   python analysis/label_map_alt.py job TISSUE DATASET SCRIPT
#   python analysis/label_map_alt.py summarise
# Output: label_map_alt/<tissue>/<dataset>/<script>.csv, label_map_alt/map_alt_all.csv, log label_map_alt/run.log
import os, sys, time, glob, runpy, subprocess
import numpy as np
import pandas as pd

# ---------------- config ----------------
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/label_map_alt'
PY = os.path.join(os.environ.get('CONDA_ENVS', os.path.expanduser('~/miniconda3/envs')), 'scib/bin/python')
JOBS = 10
THREADS_DISSECTED = 6      # kNN float32 results depend on the thread count: the saved MAP used 6 threads for the
                           # dissected datasets (extra_metrics run of 31 Oct) and extra_metrics.THREADS (4) for genes2-8
MIN_FREE_GB = 60
CHECK_TOL = 1e-9           # EBCC MAP recomputed here vs extra_metrics/extra_metrics_all.csv
SPECIES = 'orig.ident'
LABEL_COLS = {             # best annotator per species (label_sensitivity_choice.json; as the 05 label-sensitivity runs)
    'bone_marrow':   {'hsapiens': 'celltype_seurat', 'mmusculus': 'ebcc_label'},
    'fetal_liver':   {'hsapiens': 'celltype_seurat', 'mmusculus': 'ebcc_label'},
    'breast_cancer': {'hsapiens': 'celltype_seurat', 'mmusculus': 'celltype_seurat_cross'},
}
GENES = [f'genes{g}' for g in (2, 4, 6, 8)]
DATASETS = {'fetal_liver': ['dissected_mouse', *GENES], 'breast_cancer': ['dissected_mouse', *GENES],
            'bone_marrow': ['dissected', *GENES]}            # longest first
GRAPH_ONLY = ('BBKNN', 'SAMap')
# -----------------------------------------

sys.path.insert(0, f'{ROOT}/analysis')
import extra_metrics as em
THREADS_OVERRIDE = {('bone_marrow', 'dissected', '05.11_evaluate_LIGER_UINMF.py'): 4}   # failed in the 6-thread dissected
                           # run (extra_metrics/run_dissected.log) and was rescored in the 4-thread run
def threads(tissue, dataset, script):
    return THREADS_OVERRIDE.get((tissue, dataset, script), THREADS_DISSECTED if dataset.startswith('dissected') else em.THREADS)

def scripts(tissue):
    return sorted(os.path.basename(f) for f in glob.glob(f'{ROOT}/{tissue}/05.*.py')
                  if 'scIB' not in f and not any(g in os.path.basename(f) for g in GRAPH_ONLY))

def out_file(tissue, dataset, script): return f'{OUT}/{tissue}/{dataset}/{script[:-3]}.csv'

def map_two(x, labs):
    """em.mean_average_precision for several label vectors from one neighbour search"""
    n = x.shape[0]
    k = max(round(n * em.NEIGHBOR_FRAC), 1)
    ap = {name: np.empty(n) for name in labs}
    for s, nni in em._knn(x, k):
        for name, y in labs.items():
            match = y[nni[:, 1:]] == y[s:s + nni.shape[0], None]
            prec = np.cumsum(match, axis = 1) / np.arange(1, match.shape[1] + 1)
            m = match.sum(axis = 1)
            ap[name][s:s + nni.shape[0]] = np.where(m > 0, (prec * match).sum(axis = 1) / np.maximum(m, 1), 0.0)
    return {name: float(v.mean()) for name, v in ap.items()}, k

class _Captured(Exception):
    pass

def job(tissue, dataset, script):
    em.THREADS = threads(tissue, dataset, script)
    os.chdir(f'{ROOT}/{tissue}')
    sys.path.insert(0, os.getcwd())
    sys.argv = [script, dataset]
    import utils
    captured = []
    def score(adata_int, method, output_type, raw, embed = None, **kw):
        t0 = time.time()
        res = dict(tissue = tissue, dataset = dataset, script = script, method = method, output_type = output_type,
                   embed = embed, n_cells = adata_int.n_obs)
        assert output_type == 'embed', f'{script}: not an embedding method'
        assert (raw.obs_names == adata_int.obs_names).all(), 'raw and integrated cells differ'
        x = np.asarray(adata_int.obsm[embed], dtype = np.float32)
        obs = adata_int.obs if em.LABEL in adata_int.obs else raw.obs   # as analysis/extra_metrics.py
        ebcc = obs[em.LABEL].astype(str).to_numpy()
        sp = raw.obs[SPECIES].astype(str).to_numpy()                  # as utils.score builds a sensitivity label set
        alt = np.empty(len(sp), dtype = object)
        for s_, col in LABEL_COLS[tissue].items():
            alt[sp == s_] = raw.obs[col].astype(str).to_numpy()[sp == s_]
        assert set(np.unique(sp)) == set(LABEL_COLS[tissue]), 'unexpected species'
        v, k = map_two(x, {'ebcc': ebcc, 'best_annotator': alt.astype(str)})
        res.update(MAP_ebcc = v['ebcc'], MAP_best_annotator = v['best_annotator'], k_MAP = k,
                   n_labels_ebcc = len(np.unique(ebcc)), n_labels_alt = len(np.unique(alt.astype(str))),
                   metric_seconds = time.time() - t0)
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
    print(f'saved {out}: MAP ebcc {res["MAP_ebcc"]:.6f}, best annotator {res["MAP_best_annotator"]:.6f} ({res["seconds"]:.0f} s)')

def log(msg):
    t = time.strftime("%H:%M:%S")   # '=====' markers start the line, for ~/bin/watch_job.sh
    line = f'{msg} ({t})' if msg.startswith('=====') else f'{t} {msg}'
    print(line, flush = True)
    with open(f'{OUT}/run.log', 'a') as f:
        f.write(line + '\n')

def summarise():
    rows = [pd.read_csv(f) for f in glob.glob(f'{OUT}/*/*/*.csv')]
    d = pd.concat(rows)
    saved = pd.read_csv(f'{ROOT}/extra_metrics/extra_metrics_all.csv')[['tissue', 'dataset', 'method', 'MAP']]
    c = d.merge(saved, on = ['tissue', 'dataset', 'method'], how = 'left')
    c['diff'] = (c.MAP_ebcc - c.MAP).abs()
    log(f'check vs saved EBCC MAP: {c.MAP.notna().sum()} / {len(c)} matched, max |diff| {c["diff"].max():.3g}')
    d.sort_values(['tissue', 'dataset', 'method']).to_csv(f'{OUT}/map_alt_all.csv', index = False)
    bad = c[(c['diff'] > CHECK_TOL) | c.MAP.isna()]
    if len(bad):
        log(f'===== FAILED: check, {len(bad)} rows differ from the saved MAP'); print(bad.to_string()); sys.exit(1)
    log(f'wrote {OUT}/map_alt_all.csv ({len(d)} rows)')

def run():
    os.makedirs(OUT, exist_ok = True)
    jobs = [(t, d, s) for t in DATASETS for d in DATASETS[t] for s in scripts(t)]
    todo = [j for j in jobs if not os.path.exists(out_file(*j))]
    log(f'===== start: {len(todo)} jobs ({len(jobs) - len(todo)} already done)')
    tvars = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS', 'N_THREADS', 'N_CPU')
    running, failed = {}, []
    while todo or running:
        for p, j in list(running.items()):
            if p.poll() is not None:
                del running[p]
                if p.returncode != 0: failed.append(j)
                log(f'{"done  " if p.returncode == 0 else "FAILED"} {" ".join(j)}')
        while todo and len(running) < JOBS and em.free_gb() > MIN_FREE_GB:
            j = todo.pop(0)
            lf = f'{OUT}/logs/{j[0]}/{j[1]}/{j[2][:-3]}.log'
            os.makedirs(os.path.dirname(lf), exist_ok = True)
            env = dict(os.environ, **{v: str(threads(*j)) for v in tvars})
            p = subprocess.Popen([PY, os.path.abspath(__file__), 'job', *j], env = env,
                                 stdout = open(lf, 'w'), stderr = subprocess.STDOUT, cwd = ROOT)
            running[p] = j
            log(f'start  {" ".join(j)}')
        time.sleep(10)
    if failed:
        log(f'===== FAILED: {len(failed)} jobs: {failed}'); sys.exit(1)
    summarise()
    log('===== ALL DONE')

if __name__ == '__main__':
    if sys.argv[1] == 'job': job(*sys.argv[2:5])
    elif sys.argv[1] == 'run': run()
    elif sys.argv[1] == 'summarise': summarise()
