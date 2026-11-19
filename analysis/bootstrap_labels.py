# Fetal liver, dissected: the analysis/bootstrap_scores.py subsamples rescored under a second label set, plus the two
# unintegrated inputs (Unintegrated UINMF, Unintegrated scVI-latent) that analysis/bootstrap_scores.py does not score.
# Reuses bootstrap_scores.job unchanged (same 50 subsamples, paired with bootstrap/; each tissue's own 05 script
# and utils.score; no integration rerun), with the label set from LABEL_SET / LABEL_COLS as in the
# label-sensitivity rescoring. Writes only under bootstrap_labels/; bootstrap/ is read, never written.
#   python analysis/bootstrap_labels.py run      all jobs (resumable), log bootstrap_labels/run.log
#   python analysis/bootstrap_labels.py job LABELS REP SCRIPT   one job (LABELS = ebcc or best_annotator)
# Output: bootstrap_labels/<labels>/scores/<rep>/<script>.csv
import os, sys, glob, time, subprocess
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f'{ROOT}/analysis')
import bootstrap_scores as bs

TISSUE = 'fetal_liver'
OUT = f'{ROOT}/bootstrap_labels'
REPS = range(1, 51)
LABELS = {'ebcc': None,
          'best_annotator': 'hsapiens=celltype_seurat,mmusculus=ebcc_label'}   # as the 05 label-sensitivity run
JOBS = int(os.environ.get('BOOT_JOBS', '16'))
SUB_DIR = f'{ROOT}/bootstrap/{TISSUE}/subsamples'          # the original subsamples (read only)

def scripts():
    every = sorted(os.path.basename(f) for f in glob.glob(f'{ROOT}/{TISSUE}/05.*.py') if 'scIB' not in f)
    return every

def jobs():
    """best_annotator: every script; ebcc: only the scripts bootstrap/ has no scores for"""
    have = {os.path.basename(f)[:-4] + '.py' for f in glob.glob(f'{ROOT}/bootstrap/{TISSUE}/scores/1/*.csv')}
    out = [('best_annotator', r, s) for r in REPS for s in scripts()]
    out += [('ebcc', r, s) for r in REPS for s in scripts() if s not in have]
    return out

def score_path(labels, rep, script): return f'{OUT}/{labels}/scores/{rep}/{script[:-3]}.csv'

def log(msg):
    line = f'{time.strftime("%H:%M:%S")} {msg}'
    print(line, flush = True)
    with open(f'{OUT}/run.log', 'a') as f: f.write(line + '\n')

if sys.argv[1] == 'job':
    labels, rep, script = sys.argv[2], sys.argv[3], sys.argv[4]
    if LABELS[labels]:
        os.environ['LABEL_SET'], os.environ['LABEL_COLS'] = labels, LABELS[labels]
    bs.sub_file = lambda tissue, b: f'{SUB_DIR}/rep_{b:03d}.txt'
    bs.score_file = lambda tissue, r, s: score_path(labels, r, s)
    bs.job(TISSUE, rep, script)

elif sys.argv[1] == 'run':
    os.makedirs(OUT, exist_ok = True)
    todo = [j for j in jobs() if not os.path.exists(score_path(*j))]
    log(f'===== start: {len(todo)} jobs, {JOBS} at a time')
    env = dict(os.environ, **{v: str(bs.THREADS) for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                                                           'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS', 'N_THREADS', 'N_CPU')})
    running, failed = {}, []
    while todo or running:
        for p, j in list(running.items()):
            if p.poll() is not None:
                del running[p]
                if p.returncode != 0: failed.append(j)
                log(f'{"done  " if p.returncode == 0 else "FAILED"} {j[0]} {j[1]} {j[2]}')
        while todo and len(running) < JOBS and bs.free_gb() > bs.MIN_FREE_GB:
            j = todo.pop(0)
            lf = f'{OUT}/logs/{j[0]}/{j[1]}/{j[2][:-3]}.log'
            os.makedirs(os.path.dirname(lf), exist_ok = True)
            p = subprocess.Popen([bs.PY, os.path.abspath(__file__), 'job', *map(str, j)], env = env,
                                 stdout = open(lf, 'w'), stderr = subprocess.STDOUT, cwd = ROOT)
            running[p] = j
        time.sleep(10)
    if failed:
        log(f'===== FAILED: {len(failed)} jobs, e.g. {failed[:3]}')
        sys.exit(1)
    log('===== ALL DONE')
