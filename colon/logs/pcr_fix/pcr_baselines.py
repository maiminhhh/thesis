# PCR_batch baseline (variance explained by species in the unintegrated PCA) per dataset, on the 05 evaluation
# cells, with the old baseline (PCA of the per-species scale.data layer) and the fixed one (joint_pca).
# PCR_batch = max(0, (before - after) / before) and `after` depends only on the method, so the old score
# of every method follows from the new 05 results: see compare_pcr.py.
import os, sys, warnings
warnings.filterwarnings('ignore')
COLON = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))   # the colon folder
sys.path.insert(0, COLON); os.chdir(COLON)
import scanpy as sc, pandas as pd
from scib.metrics.pcr import pcr
from utils import joint_pca, eval_cells, BATCH
from integration_datasets import DATASETS
rows = []
for d, ds in DATASETS.items():
    raw = sc.read_h5ad(f"{ds['concat']}.h5ad")
    raw = raw[raw.obs_names.intersection(eval_cells(ds['out']), sort=False)].copy()
    old = raw.copy()
    data_layer = next((l for l in old.layers if l.startswith('data')), None)
    scale_layer = next((l for l in old.layers if l.startswith('scale.data')), None)
    old.X = old.layers[data_layer]
    sc.pp.pca(old, layer=scale_layer)
    joint_pca(raw)
    b_old = pcr(old, covariate=BATCH, linreg_method='numpy')
    b_new = pcr(raw, covariate=BATCH, linreg_method='numpy')
    rows.append(dict(dataset=d, n_cells=raw.n_obs, n_genes=raw.n_vars, before_old=b_old, before_new=b_new))
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv('logs/pcr_fix/pcr_baselines.csv', index=False)
