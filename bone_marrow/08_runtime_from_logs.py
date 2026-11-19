# Approximate runtime of every 04 method, taken from the run logs (no new runs):
#   python 08_runtime_from_logs.py   ->  data/08_timing/runtime_from_logs.csv
# Wall-clock job time on a shared server, INCLUDING reading inputs, UMAP and saving
# outputs, with up to 6 CPU jobs and 2 GPU jobs running at once. Good for tiers
# (minutes vs tens of minutes), not for ranking methods within a tier.
#   full                : logs/pipeline.log (last complete run; the full data did not change)
#   dissected, genes2-8 : logs/rerun_dissect.log (current dissection)
#   SAMap dissected     : logs/rerun_samap_dissected.log (rerun with the refreshed input)
#   scVI latent (04.18) : logs/scvi_latent/gpu.log (genes8: its own log's file times)
#   04.19-04.24         : file creation -> last write time of each job log (no timestamps logged)
# minutes_excl_umap_init removes uwot's spectral initialisation, which stalled for up to
# ~4 h on some datasets (plotting only, not integration).
import os, re, glob
import datetime as dt
import pandas as pd
import anndata as ad
from datasets import DATASETS

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
# the one-off step each method needs before it can run (its own time is in its own row)
NEEDS = {'04.08': '04.07 scVI', **{k: '04.11 LIGER UINMF (H matrix)' for k in ['04.12', '04.13', '04.14', '04.15', '04.16', '04.17']},
         **{k: '04.18 scVI latent' for k in ['04.19', '04.20', '04.21', '04.22', '04.23', '04.24']}}

def runner_times(log, datasets):
    """last start/done pair per job in a run_03_04.sh-style log -> {(script, dataset): minutes}"""
    start, out = {}, {}
    for line in open(log):
        m = re.match(r'(\d\d:\d\d:\d\d) (start|done ) (04\.\d\d)_(\w+)', line)
        if not m or m.group(4) not in datasets:
            continue
        t, key = dt.datetime.strptime(m.group(1), '%H:%M:%S'), (m.group(3), m.group(4))
        if m.group(2) == 'start':
            start[key] = t
        elif key in start:
            out[key] = ((t - start.pop(key)).total_seconds() % 86400) / 60   # jobs < 24 h
    return out

rows = []
def add(times, source):
    for (script, s), minutes in times.items():
        rows.append(dict(script = script, dataset = s, minutes = round(minutes, 2), source = source))

add(runner_times('logs/pipeline.log', {'full'}), 'logs/pipeline.log')
add(runner_times('logs/rerun_dissect.log', {'dissected', 'genes2', 'genes4', 'genes6', 'genes8'}),
    'logs/rerun_dissect.log')
# SAMap dissected: rerun after its per-species input was refreshed (supersedes rerun_dissect)
steps = [(dt.datetime.strptime(m.group(2), '%b %d %H:%M:%S %Z %Y'), m.group(1))
         for m in (re.match(r'===== (.+?) \w{3} (\w{3} +\d+ [\d:]+ \w+ \d{4})', l)
                   for l in open('logs/rerun_samap_dissected.log')) if m]
for (t0, name), (t1, _) in zip(steps, steps[1:]):
    if name.startswith('04.02 SAMap dissected'):
        rows = [r for r in rows if (r['script'], r['dataset']) != ('04.02', 'dissected')]
        add({('04.02', 'dissected'): (t1 - t0).total_seconds() / 60}, 'logs/rerun_samap_dissected.log')
# 04.18 on the GPUs
start = {}
for line in open('logs/scvi_latent/gpu.log'):
    m = re.match(r'(\d\d:\d\d:\d\d) (start|done ) 04\.18 (\w+)', line)
    if m:
        t = dt.datetime.strptime(m.group(1), '%H:%M:%S')
        if m.group(2) == 'start': start[m.group(3)] = t
        else: add({('04.18', m.group(3)): ((t - start[m.group(3)]).total_seconds() % 86400) / 60},
                  'logs/scvi_latent/gpu.log')
# job logs without timestamps: file creation -> last write. genes8 04.19/20/21/23 were rerun
# into the same file after a failed first run, so their creation time is wrong: skipped
for f in glob.glob('logs/scvi_latent/04.*_*.log'):
    script, s = re.match(r'.*/(04\.\d\d)_(\w+)\.log', f).groups()
    if script == '04.18' and s != 'genes8':
        continue
    if s == 'genes8' and script in ('04.19', '04.20', '04.21', '04.23'):
        continue
    st = os.stat(f)
    born = getattr(st, 'st_birthtime', None) or int(os.popen(f'stat -c %W "{f}"').read())
    if born > 0:
        add({(script, s): (st.st_mtime - born) / 60}, 'log file times')

def umap_init_minutes(script, s):
    """time uwot spent on its spectral initialisation (R jobs; it can stall for hours on
    some datasets, e.g. genes2): from 'Initializing from normalized Laplacian' to the next
    timestamped line. The UMAP is only for plotting, so this is not integration time."""
    for f in (f'logs/scvi_latent/{script}_{s}.log', f'logs/03_04/{script}_{s}.log'):
        if not os.path.exists(f):
            continue
        stamps = [(dt.datetime.strptime(m.group(1), '%H:%M:%S'), line)
                  for line in open(f, errors = 'ignore')
                  for m in [re.match(r'(\d\d:\d\d:\d\d) ', line)] if m]
        for (t0, line), (t1, _) in zip(stamps, stamps[1:]):
            if 'Initializing from normalized Laplacian' in line:
                return ((t1 - t0).total_seconds() % 86400) / 60
        return 0.0
    return float('nan')

df = pd.DataFrame(rows)
df['umap_init_minutes'] = [round(umap_init_minutes(k, s), 2) for k, s in zip(df.script, df.dataset)]
df['minutes_excl_umap_init'] = (df.minutes - df.umap_init_minutes.fillna(0)).round(2)
df['method'] = df.script.map(lambda k: METHODS[k][0])
df['input'] = df.script.map(lambda k: METHODS[k][1])
df['device'] = df.script.map(lambda k: 'gpu' if k in GPU else 'cpu')
df['needs'] = df.script.map(NEEDS).fillna('')
n_cells = {s: ad.read_h5ad(d['concat'], backed = 'r').n_obs for s, d in DATASETS.items()}
df['n_cells'] = df.dataset.map(n_cells)
df['note'] = 'wall-clock incl. I/O, UMAP, saving; shared server, up to 6 CPU + 2 GPU jobs in parallel'
df.loc[df.script.isin(['04.19', '04.20', '04.21', '04.22', '04.23', '04.24']), 'note'] += \
    '; ran in a separate, busier batch (from log file times)'
order = ['full', 'dissected', 'genes2', 'genes4', 'genes6', 'genes8']
df = df.sort_values(['script', 'dataset'], key = lambda c: c.map({s: i for i, s in enumerate(order)}) if c.name == 'dataset' else c)
cols = ['script', 'method', 'input', 'device', 'needs', 'dataset', 'n_cells', 'minutes',
        'umap_init_minutes', 'minutes_excl_umap_init', 'source', 'note']
os.makedirs('data/08_timing', exist_ok = True)
df[cols].to_csv('data/08_timing/runtime_from_logs.csv', index = False)
print(df.pivot(index = 'method', columns = 'dataset', values = 'minutes_excl_umap_init').reindex(columns = order).round(1).to_string())
print(f'{len(df)} rows -> data/08_timing/runtime_from_logs.csv')
