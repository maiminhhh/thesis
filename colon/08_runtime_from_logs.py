# Approximate runtime of every 04 method, taken from the run logs (no new runs; ../EBCC_PIPELINE.md 8b):
#   python 08_runtime_from_logs.py   ->  data/08_timing/runtime_from_logs.csv
# Wall-clock job time on a shared server, INCLUDING reading inputs, UMAP and saving outputs, with
# up to 4-6 CPU jobs and 2 GPU jobs running at once. Good for tiers (minutes vs tens of minutes),
# not for ranking methods within a tier. Same columns as bone_marrow's, datasets named full /
# dissected / genes2-8, so the tissues can be plotted together.
# Sources (latest complete run on the current data; the 03 rerun of 2026-10-24 changed labels only,
# cells and genes identical, so the 04 outputs of the run before it are current):
#   logs/04/status.txt                 : run_04.sh, 2026-10-23/24 (minutes per job, whole minutes)
#   logs/04/status_samap_rerun.txt     : SAMap dissected/genes4/genes8 (rerun after GPU out-of-memory)
#   logs/04/status_fewgenes_rerun.txt  : LIGER genes6/genes8, RPCA UINMF genes8 (rerun with capped k / dims)
#   logs/04/status_scvi.txt            : 04.18 scVI latent (and 04.19-04.24 once run)
#   scANVI (04.08)                     : rerun 2026-10-24 on the final labels, two GPU lanes one job after
#                                        another from the "03 rerun done" line; its logs overwrote the first
#                                        run's (creation time is from that run), so each job's time is the
#                                        gap between consecutive log writes in its lane
# minutes_excl_umap_init removes uwot's spectral initialisation ("Initializing from normalized
# Laplacian"; only 04.11 uses it here, the rest use a PCA init): plotting only, not integration.
import os, re, datetime as dt
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
NAME = {'full_mouse': 'full', 'dissected_mouse': 'dissected', **{f'genes{g}': f'genes{g}' for g in (2, 4, 6, 8)}}

rows = {}
def add(script, ds, minutes, source):   # later sources replace earlier ones (reruns on current inputs)
    rows[(script, ds)] = dict(script = script, dataset = NAME[ds], minutes = round(minutes, 2), source = source)

def status_file(path):
    """run_job lines: '<script> <dataset> exit=N  M min'; successful jobs only"""
    if not os.path.exists(path): return
    for line in open(path):
        m = re.match(r'(04\.\d\d)\S*\s+(\S+)\s+exit=(\d+)\s+(\d+) min', line)
        if m and m.group(3) == '0' and m.group(2) in NAME:
            add(m.group(1), m.group(2), float(m.group(4)), path)

status_file('logs/04/status.txt')
status_file('logs/04/status_samap_rerun.txt')
status_file('logs/04/status_fewgenes_rerun.txt')
status_file('logs/04/status_scvi.txt')

# scANVI rerun: two GPU lanes from "03 rerun done", jobs back to back -> gaps between log write times
t0 = None
for line in open('logs/pipeline_status.txt'):
    if '03 rerun done' in line:
        t0 = dt.datetime.strptime(line[:19], '%Y-%m-%d %H:%M:%S').timestamp()
for lane in (['full_mouse', 'genes2', 'genes6'], ['dissected_mouse', 'genes4', 'genes8']):
    prev = t0
    for ds in lane:
        end = os.path.getmtime(f'logs/04/04.08_integrate_scANVI__{ds}.log')
        add('04.08', ds, (end - prev) / 60, 'scANVI rerun lane: gap between job-log writes')
        prev = end

def umap_init_minutes(script, ds):
    """uwot spectral initialisation: from 'Initializing from normalized Laplacian' to the next timestamped line"""
    f = next((f for f in [f'logs/04/{s}' for s in os.listdir('logs/04')] if re.match(rf'logs/04/{re.escape(script)}_.*__{ds}\.log$', f)), None)
    if f is None: return 0.0
    stamps = [(dt.datetime.strptime(m.group(1), '%H:%M:%S'), line) for line in open(f, errors = 'ignore')
              for m in [re.match(r'(\d\d:\d\d:\d\d) ', line)] if m]
    for (a, line), (b, _) in zip(stamps, stamps[1:]):
        if 'Initializing from normalized Laplacian' in line:
            return ((b - a).total_seconds() % 86400) / 60
    return 0.0

df = pd.DataFrame(rows.values())
inv = {v: k for k, v in NAME.items()}
df['umap_init_minutes'] = [round(umap_init_minutes(k, inv[s]), 2) for k, s in zip(df.script, df.dataset)]
df['minutes_excl_umap_init'] = (df.minutes - df.umap_init_minutes).round(2)
df['method'] = df.script.map(lambda k: METHODS[k][0])
df['input'] = df.script.map(lambda k: METHODS[k][1])
df['device'] = df.script.map(lambda k: 'gpu' if k in GPU else 'cpu')
df['needs'] = df.script.map(NEEDS).fillna('')
n_cells = {NAME[s]: ad.read_h5ad(f"{d['concat']}.h5ad", backed = 'r').n_obs for s, d in DATASETS.items()}
df['n_cells'] = df.dataset.map(n_cells)
df['note'] = 'wall-clock incl. I/O, UMAP, saving; shared server, up to 4-6 CPU + 2 GPU jobs in parallel'
df.loc[df.source.str.contains('status.txt$') & (df.device == 'cpu'), 'note'] += '; whole minutes'
# outliers: far slower than the same method on the other datasets, not explained by cell numbers
# (full has 1.6x the cells of dissected). Kept as logged (no clean way to correct them from the logs), flagged here.
OUTLIERS = {('04.07', 'dissected'): 'OUTLIER: ran under heavy load, 101 s/epoch vs 38 s/epoch on genes2 (same cells)',
            ('04.07', 'full'):      'OUTLIER: ran under heavy load, 194 s/epoch',
            ('04.02', 'full'):      'OUTLIER: 4-6x the other datasets (shared GPU / busy server)',
            ('04.05', 'full'):      'OUTLIER: ~4-5x the other datasets',
            ('04.15', 'full'):      'OUTLIER: ~5-8x the other datasets'}
for (k, s), msg in OUTLIERS.items():
    df.loc[(df.script == k) & (df.dataset == s), 'note'] += '; ' + msg
order = ['full', 'dissected', 'genes2', 'genes4', 'genes6', 'genes8']
df = df.sort_values(['script', 'dataset'], key = lambda c: c.map({s: i for i, s in enumerate(order)}) if c.name == 'dataset' else c)
cols = ['script', 'method', 'input', 'device', 'needs', 'dataset', 'n_cells', 'minutes',
        'umap_init_minutes', 'minutes_excl_umap_init', 'source', 'note']
os.makedirs('data/08_timing', exist_ok = True)
df[cols].to_csv('data/08_timing/runtime_from_logs.csv', index = False)
pd.set_option('display.width', 200)
print(df.pivot(index = 'method', columns = 'dataset', values = 'minutes_excl_umap_init').reindex(columns = order).round(0).to_string())
print(df[df.umap_init_minutes > 0][['method', 'dataset', 'minutes', 'umap_init_minutes']].to_string(index = False))
print(f'{len(df)} rows -> data/08_timing/runtime_from_logs.csv')
