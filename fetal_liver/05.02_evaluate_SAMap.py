print('1. Evaluate SAMap - Setup')
import os
N_CPU = 32
os.environ.setdefault("OMP_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
os.environ.setdefault("MKL_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
os.environ.setdefault("NUMEXPR_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
import scanpy as sc
from utils import score, align, eval_cells
from integration_datasets import DATASETS

method = 'SAMap'
for g, d in DATASETS.items():
    print(f'2. Evaluate {method} - {g}')
    raw = sc.read_h5ad(f"{d['concat']}.h5ad")
    adata = sc.read_h5ad(f"{d['out']}04.02_SAMap.h5ad")
    species = ['mmusculus' if s == 'mm' else 'hsapiens' for s in adata.obs['orig.ident'].astype(str)]
    adata.obs_names = [f'{s}_{c}' for s, c in zip(species, adata.obs_names)]
    # SAMap's output lacks a few evaluation cells (<= 9 mouse cells): SAMap alone is scored
    # on the evaluation cells it has
    cells = eval_cells(d['out'])
    print(f'SAMap: {len(cells.difference(adata.obs_names))} evaluation cells not in the SAMap output, left out')
    adata, raw = align(adata, raw, cells.intersection(adata.obs_names, sort = False))
    sc.pp.pca(adata)
    scib_scores = score(adata_int = adata,
                        method = method,
                        output_type = 'knn',
                        raw = raw)

    print(f'3. Evaluate {method} - Save score - {g}')
    scib_scores.to_csv(f"{d['scores']}05.02_SAMap.csv")
    print(f'Evaluate {method} - {g} - Done!')
