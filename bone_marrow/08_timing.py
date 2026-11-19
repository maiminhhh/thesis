# Runtime of one Python integration method, integration step only (no UMAP, no saving):
#   python 08_timing.py <method> <dataset> <repeat>
# Inputs are read before the clock starts. PCA-based methods include their PCA (as the R
# helpers do). UINMF / scVI-latent variants read their input from disk; making it is timed
# separately ('LIGER UINMF' in 08_timing.R, 'scVI latent' here). scANVI starts from the
# saved 04.07 scVI model, so its full cost is scVI + scANVI.
# Appends one row to data/08_timing/timing.csv (written by run_timing.sh).
import os
N = int(os.environ.get('N_THREADS', 8))
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[v] = str(N)
import sys, time, gc
import numpy as np
import pandas as pd
import scanpy as sc
from datasets import DATASETS

method, name, rep = sys.argv[1], sys.argv[2], int(sys.argv[3])
d = DATASETS[name]
GPU = {'scVI', 'scANVI', 'scVI latent', 'SAMap'}

def read_basis(adata, variant):
    """X_basis for the UINMF / scVI-latent variants, cells as in 04.15/17 and 04.22/24"""
    if variant == 'UINMF':
        H = pd.read_csv(f"{d['out']}uinmf_H_raw.csv", index_col = 0)
        usable = H.index[np.isfinite(H).all(axis = 1) & (H.abs().sum(axis = 1) > 0)]
        adata = adata[adata.obs_names.intersection(usable, sort = False)].copy()
    else:
        H = pd.read_csv(f"{d['out']}scvi_latent.csv", index_col = 0)
    adata.obsm['X_basis'] = H.loc[adata.obs_names].to_numpy(dtype = float)
    return adata

base, _, variant = method.partition(' ')
if method == 'SAMap':
    from samap import SAMAP
    from samap.utils import prepend_var_prefix
    from samap.sam import SAM
    inputs = {}
    for sp in ('human', 'mouse'):
        a = sc.read_h5ad(d[sp]); a.var_names = a.var['names']
        for key in ['data', 'scale.data', 'X_disp', 'X_counts']:
            a.layers.pop(key, None)
        inputs[sp] = a
    n_cells = sum(a.n_obs for a in inputs.values())
else:
    adata = sc.read_h5ad(d['concat'])
    n_cells = adata.n_obs
    if base in ('Scanorama', 'BBKNN') and variant:
        adata = read_basis(adata, variant)
    if method in ('scVI', 'scANVI', 'scVI latent'):
        import scvi, torch
        torch.set_num_threads(N)
        scvi.settings.seed = 1
        adata.layers['counts'] = adata.X
    if base == 'BBKNN':
        import sklearn.metrics, sklearn.neighbors
        sklearn.neighbors.DistanceMetric = sklearn.metrics.DistanceMetric
        sys.modules['sklearn.neighbors.DistanceMetric'] = sklearn.metrics.DistanceMetric
        import bbknn
gc.collect()

start = time.perf_counter()
if method == 'Scanorama':
    sc.pp.pca(adata, layer = 'scale.data_Xenium')
    sc.external.pp.scanorama_integrate(adata, key = 'orig.ident', basis = 'X_pca', adjusted_basis = 'X_out')
elif base == 'Scanorama':
    sc.external.pp.scanorama_integrate(adata, key = 'orig.ident', basis = 'X_basis', adjusted_basis = 'X_out')
elif method == 'BBKNN':
    sc.pp.pca(adata, layer = 'scale.data_Xenium')
    bbknn.bbknn(adata, batch_key = 'orig.ident')
elif base == 'BBKNN':
    bbknn.bbknn(adata, batch_key = 'orig.ident', use_rep = 'X_basis')
elif method == 'scVI':
    scvi.model.SCVI.setup_anndata(adata, layer = 'counts', batch_key = 'orig.ident')
    m = scvi.model.SCVI(adata); m.train(); m.get_latent_representation()
elif method == 'scVI latent':
    scvi.model.SCVI.setup_anndata(adata, layer = 'counts')
    m = scvi.model.SCVI(adata, n_latent = 20); m.train(); m.get_latent_representation()
elif method == 'scANVI':
    m = scvi.model.SCVI.load(f"{d['out']}04.07_scVI_model", adata = adata)
    m = scvi.model.SCANVI.from_scvi_model(m, labels_key = 'ebcc_label', unlabeled_category = 'Unknown')
    m.train(); m.get_latent_representation()
elif method == 'SAMap':
    sams = {}
    for sp, prefix in (('human', 'hs'), ('mouse', 'mm')):
        s = SAM(inputs[sp]); prepend_var_prefix(s, prefix)
        s.preprocess_data(); s.run(); sams[prefix] = s
    sm = SAMAP(sams, f_maps = 'data/04_integration_resources/maps/', backend = 'cuda', save_processed = False)
    sm.run()
else:
    raise SystemExit(f'unknown method: {method}')
sec = time.perf_counter() - start

row = pd.DataFrame([dict(method = method, dataset = name, rep = rep, seconds = round(sec, 2),
                         device = 'gpu' if method in GPU else 'cpu', threads = N, n_cells = n_cells,
                         load1 = float(open('/proc/loadavg').read().split()[0]),
                         time = time.strftime('%Y-%m-%d %H:%M:%S'))])
f = 'data/08_timing/timing.csv'
os.makedirs(os.path.dirname(f), exist_ok = True)
row.to_csv(f, mode = 'a', header = not os.path.exists(f), index = False)
print(row.to_string(index = False))
