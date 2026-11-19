# Approximate runtime of every 04 method, taken from the run logs (no new runs):
#   python 08_runtime_from_logs.py   ->  data/08_timing/runtime_from_logs.csv
# Wall-clock job time on a shared server, INCLUDING reading inputs, UMAP and saving
# outputs, with several CPU jobs (4 threads each) and 2 GPU jobs running at once. Good for
# tiers (minutes vs tens of minutes), not for ranking methods within a tier.
# Sources (fetal_liver; same columns as bone_marrow/08_runtime_from_logs.py):
#   logs/04/status.txt        job minutes per run (run_04.sh / run_04_scvi.sh / the v2 driver).
#                             genes2-8: only runs after the renormalised 03 (2026-10-23 14:46),
#                             except SAMap (raw counts only, unchanged by the renormalisation);
#                             full / dissected: their 03 data did not change, so earlier runs count.
#                             Later lines supersede earlier ones.
#   logs/04_uinmf/status.txt  04.12-04.17 on full / dissected (2026-10-22)
#   job log file times        (creation -> last write) where no minutes were logged: 04.11 on
#                             full / dissected (logs/04_rerun/, 2026-10-22) and SAMap dissected
# minutes_excl_umap_init removes uwot's spectral initialisation ("Initializing from normalized
# Laplacian"; plotting only, not integration), read from the timestamps in each job log.
import os, re
import datetime as dt
import pandas as pd
import anndata as ad
from integration_datasets import DATASETS

METHODS = {
    '04.01': ('FastMNN', 'PCA'), '04.02': ('SAMap', 'own'), '04.03': ('CCA', 'PCA'),
    '04.04': ('Harmony', 'PCA'), '04.05': ('Scanorama', 'PCA'), '04.06': ('RPCA', 'PCA'),
    '04.07': ('scVI', 'own'), '04.08': ('scANVI', 'own'), '04.09': ('BBKNN', 'PCA'),
    '04.10': ('LIGER', 'own'), '04.11': ('LIGER UINMF', 'UINMF'),
    '04.12': ('FastMNN UINMF', 'UINMF'), '04.13': ('CCA UINMF', 'UINMF'),
    '04.14': ('Harmony UINMF', 'UINMF'), '04.15': ('Scanorama UINMF', 'UINMF'),
    '04.16': ('RPCA UINMF', 'UINMF'), '04.17': ('BBKNN UINMF', 'UINMF'),
    '04.18': ('scVI latent (input)', 'scVI-latent'),
    '04.19': ('FastMNN scVI-latent', 'scVI-latent'), '04.20': ('CCA scVI-latent', 'scVI-latent'),
    '04.21': ('Harmony scVI-latent', 'scVI-latent'), '04.22': ('Scanorama scVI-latent', 'scVI-latent'),
    '04.23': ('RPCA scVI-latent', 'scVI-latent'), '04.24': ('BBKNN scVI-latent', 'scVI-latent')}
GPU = {'04.02', '04.07', '04.08', '04.18'}
NEEDS = {'04.08': '04.07 scVI', **{k: '04.11 LIGER UINMF (H matrix)' for k in ['04.12', '04.13', '04.14', '04.15', '04.16', '04.17']},
         **{k: '04.18 scVI latent' for k in ['04.19', '04.20', '04.21', '04.22', '04.23', '04.24']}}
# dataset names as in bone_marrow's table
NAME = {'full_mouse': 'full', 'dissected_mouse': 'dissected', 'genes2': 'genes2', 'genes4': 'genes4',
        'genes6': 'genes6', 'genes8': 'genes8'}
CURRENT_03 = dt.datetime(2026, 10, 23, 14, 46, 43)   # first 04 run on the renormalised gene-drop data

def status_times(path, log_dir):
    """{(script, dataset): (minutes, job log)} from a status file; the last exit=0 line wins.
    genes2-8 lines count only after the CURRENT_03 run header"""
    out, run_time = {}, None
    for line in open(path):
        h = re.match(r'# run (\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)', line)
        if h:
            run_time = dt.datetime.strptime(h.group(1), '%Y-%m-%d %H:%M:%S')
            continue
        m = re.match(r'(04\.\d\d)_(\S+)\s+(\S+)\s+exit=0\s+(\S+) min', line)
        if not m:
            continue
        script, d, minutes = m.group(1), m.group(3), m.group(4)
        # SAMap reads raw counts only, which the renormalisation did not change: its earlier runs count
        if d.startswith('genes') and script != '04.02' and (run_time is None or run_time < CURRENT_03):
            continue
        log = os.path.join(log_dir, f'{m.group(1)}_{m.group(2)[:-len(os.path.splitext(m.group(2))[1])]}__{d}.log')
        out[(script, d)] = (None if minutes == '?' else float(minutes), log)
    return out

def file_minutes(f):
    st = os.stat(f)
    born = int(os.popen(f'stat -c %W "{f}"').read())
    return (st.st_mtime - born) / 60 if born > 0 else float('nan')

rows = []
def add(script, d, minutes, source, log):
    rows.append(dict(script = script, dataset = d, minutes = round(minutes, 2), source = source, log = log))

times = status_times('logs/04_uinmf/status.txt', 'logs/04_uinmf')
times = {k: v for k, v in times.items() if k[1] in ('full_mouse', 'dissected_mouse')}
times.update(status_times('logs/04/status.txt', 'logs/04'))
for (script, d), (minutes, log) in times.items():
    if minutes is None:     # logged without minutes: the job log's own file times
        add(script, d, file_minutes(log), 'log file times', log)
    else:
        add(script, d, minutes, os.path.dirname(log) + '/status.txt', log)
# LIGER UINMF on full / dissected: the run whose H matrix the UINMF variants used (logs/04_rerun)
for d in ('full_mouse', 'dissected_mouse'):
    log = f'logs/04_rerun/04.11_integrate_LIGER_UINMF__{d}.log'
    rows = [r for r in rows if (r['script'], r['dataset']) != ('04.11', d)]
    add('04.11', d, file_minutes(log), 'log file times', log)

def umap_init_minutes(log):
    """uwot's spectral initialisation: from 'Initializing from normalized Laplacian' to the next
    timestamped line (0 if the job used another initialisation)"""
    if not os.path.exists(log):
        return float('nan')
    stamps = [(dt.datetime.strptime(m.group(1), '%H:%M:%S'), line)
              for line in open(log, errors = 'ignore')
              for m in [re.match(r'(\d\d:\d\d:\d\d) ', line)] if m]
    for (t0, line), (t1, _) in zip(stamps, stamps[1:]):
        if 'Initializing from normalized Laplacian' in line:
            return ((t1 - t0).total_seconds() % 86400) / 60
    return 0.0

df = pd.DataFrame(rows)
df['umap_init_minutes'] = [round(umap_init_minutes(f), 2) for f in df.log]
df['minutes_excl_umap_init'] = (df.minutes - df.umap_init_minutes.fillna(0)).round(2)
df['method'] = df.script.map(lambda k: METHODS[k][0])
df['input'] = df.script.map(lambda k: METHODS[k][1])
df['device'] = df.script.map(lambda k: 'gpu' if k in GPU else 'cpu')
df['needs'] = df.script.map(NEEDS).fillna('')
n_cells = {d: ad.read_h5ad(f"{v['concat']}.h5ad", backed = 'r').n_obs for d, v in DATASETS.items()}
df['n_cells'] = df.dataset.map(n_cells)
df['note'] = 'wall-clock incl. I/O, UMAP, saving; shared server, several CPU jobs (4 threads) + 2 GPU jobs in parallel'
df.loc[(df.dataset == 'genes8') & df.script.isin(['04.19', '04.20', '04.21', '04.22', '04.23', '04.24']), 'note'] += \
    '; smoke test, 6 jobs at once (about twice the other datasets)'
df.loc[df.source == 'log file times', 'note'] += '; from the job log file times'
df['dataset'] = df.dataset.map(NAME)
order = list(NAME.values())
df = df.sort_values(['script', 'dataset'], key = lambda c: c.map({s: i for i, s in enumerate(order)}) if c.name == 'dataset' else c)
cols = ['script', 'method', 'input', 'device', 'needs', 'dataset', 'n_cells', 'minutes',
        'umap_init_minutes', 'minutes_excl_umap_init', 'source', 'note']
os.makedirs('data/08_timing', exist_ok = True)
df[cols].to_csv('data/08_timing/runtime_from_logs.csv', index = False)
print(df.pivot(index = 'method', columns = 'dataset', values = 'minutes_excl_umap_init')
        .reindex(index = [METHODS[k][0] for k in sorted(METHODS)], columns = order).round(1).to_string())
print(df[df.umap_init_minutes > 0][['script', 'dataset', 'minutes', 'umap_init_minutes']].to_string())
print(f'{len(df)} rows -> data/08_timing/runtime_from_logs.csv')
