print('1. Evaluate SAMap - Setup')
import os
N_CPU = int(os.environ.get('N_THREADS', 32))
os.environ["OMP_NUM_THREADS"] = str(N_CPU)
os.environ["MKL_NUM_THREADS"] = str(N_CPU)
os.environ["OPENBLAS_NUM_THREADS"] = str(N_CPU)
os.environ["NUMEXPR_NUM_THREADS"] = str(N_CPU)
import scanpy as sc
from utils import score, align, eval_cells
from integration_datasets import DATASETS

method = 'SAMap'
for d, ds in DATASETS.items():
    print(f'2. Evaluate {method} - {d}')
    raw = sc.read_h5ad(f"{ds['concat']}.h5ad")
    adata = sc.read_h5ad(f"{ds['out']}04.02_SAMap.h5ad")
    species = ['mmusculus' if s == 'mm' else 'hsapiens' for s in adata.obs['orig.ident'].astype(str)]
    adata.obs_names = [f'{s}_{c}' for s, c in zip(species, adata.obs_names)]
    adata, raw = align(adata, raw, eval_cells(ds["out"]))
    sc.pp.pca(adata)
    scib_scores = score(adata_int = adata,
                        method = method,
                        output_type = 'knn',
                        raw = raw)

    print(f'3. Evaluate {method} - Save score - {d}')
    scib_scores.to_csv(f"{ds['scores']}05.02_SAMap.csv")
    print(f'Evaluate {method} - {d} - Done!')
