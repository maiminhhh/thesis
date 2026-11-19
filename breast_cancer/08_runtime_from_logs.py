# Approximate runtime of every 04 method, taken from the run logs (no new runs; ../EBCC_PIPELINE.md 8b):
#   python 08_runtime_from_logs.py   ->  data/08_timing/runtime_from_logs.csv
# Wall-clock job time on a shared server, INCLUDING reading inputs, UMAP and saving
# outputs, with several CPU and GPU jobs running at once. Good for tiers
# (minutes vs tens of minutes), not for ranking methods within a tier.
# Adapted from bone_marrow/08_runtime_from_logs.py. This tissue's runners (run_04.sh, pipeline/run_04scvi.sh,
# pipeline/smoke_scvi_latent.sh) log no start times, so every job is timed from its own log file:
# creation (stat %W) -> last write.
#   04.01-04.17, all datasets : logs/04/<script>_<dataset>.log (run_04.sh, 2026-10-24/25, on the current 03 data)
#   04.18-04.24               : logs/04scvi/<step>_<dataset>.log (genes8: smoke test; others: run_04scvi.sh)
# Skipped: the six jobs repaired by a second run into the same log file (its creation time is from the failed first
# run): 04.11 and 04.02 full_mouse, 04.01/04.03/04.04/04.06 genes6. Their failed first attempts were renamed
# (*.umap_hang.log, *.umap_crash.log, *.gpu_oom.log) and are not used either.
# umap_init_minutes: the R job logs have no timestamps, so uwot's spectral initialisation cannot be timed (NaN, and
# minutes_excl_umap_init = minutes). genes6 R jobs use init = 'pca' (umap_helper.R), so they have no such step (0);
# Python jobs use scanpy's UMAP, not uwot (0).
import os, re, glob
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
# repaired by a second run into the same log file: creation time is from the first run
OVERWRITTEN = {('04.11', 'full_mouse'), ('04.02', 'full_mouse'),
               ('04.01', 'genes6'), ('04.03', 'genes6'), ('04.04', 'genes6'), ('04.06', 'genes6')}
R_SCRIPTS = {'04.01', '04.03', '04.04', '04.06', '04.10', '04.11', '04.12', '04.13', '04.14', '04.16',
             '04.19', '04.20', '04.21', '04.23'}

def file_minutes(f):
    """log file creation -> last write, in minutes (None if the filesystem has no creation time)"""
    born = int(os.popen(f'stat -c %W "{f}"').read())
    return (os.stat(f).st_mtime - born) / 60 if born > 0 else None

rows, skipped = [], []
logs = [(f, 'logs/04') for f in glob.glob('logs/04/04.*_*.log')] + [(f, 'logs/04scvi') for f in glob.glob('logs/04scvi/04.*_*.log')]
for f, src in logs:
    m = re.match(r'.*/(04\.\d\d)_(?:[A-Za-z_]+_)?(full_mouse|dissected_mouse|genes[2468])\.log$', f)
    if not m:      # renamed failed attempts (*.umap_hang.log etc.)
        continue
    script, s = m.groups()
    if (script, s) in OVERWRITTEN:
        skipped.append(f)
        continue
    minutes = file_minutes(f)
    if minutes is not None:
        rows.append(dict(script = script, dataset = s, minutes = round(minutes, 2), source = f'{src} (log file times)'))

df = pd.DataFrame(rows)
df['umap_init_minutes'] = [0.0 if (k not in R_SCRIPTS or s == 'genes6') else float('nan')
                           for k, s in zip(df.script, df.dataset)]
df['minutes_excl_umap_init'] = (df.minutes - df.umap_init_minutes.fillna(0)).round(2)
df['method'] = df.script.map(lambda k: METHODS[k][0])
df['input'] = df.script.map(lambda k: METHODS[k][1])
df['device'] = df.script.map(lambda k: 'gpu' if k in GPU else 'cpu')
df['needs'] = df.script.map(NEEDS).fillna('')
n_cells = {s: ad.read_h5ad(f'data/{s}/03_concatenated_data.h5ad', backed = 'r').n_obs for s in DATASETS}
df['n_cells'] = df.dataset.map(n_cells)
df['note'] = 'wall-clock incl. I/O, UMAP, saving; shared server, several CPU + GPU jobs in parallel'
df.loc[df.umap_init_minutes.isna(), 'note'] += '; uwot spectral init not timed (no timestamps in the log)'
df.loc[df.script.isin(['04.18', '04.19', '04.20', '04.21', '04.22', '04.23', '04.24']), 'note'] += \
    '; ran in a separate batch (2026-10-26)'
order = ['full_mouse', 'dissected_mouse', 'genes2', 'genes4', 'genes6', 'genes8']
df = df.sort_values(['script', 'dataset'], key = lambda c: c.map({s: i for i, s in enumerate(order)}) if c.name == 'dataset' else c)
cols = ['script', 'method', 'input', 'device', 'needs', 'dataset', 'n_cells', 'minutes',
        'umap_init_minutes', 'minutes_excl_umap_init', 'source', 'note']
os.makedirs('data/08_timing', exist_ok = True)
df[cols].to_csv('data/08_timing/runtime_from_logs.csv', index = False)
print(df.pivot(index = 'method', columns = 'dataset', values = 'minutes_excl_umap_init').reindex(columns = order).round(1).to_string())
print(f'{len(df)} rows -> data/08_timing/runtime_from_logs.csv; skipped (overwritten by a repair run): {len(skipped)}')
