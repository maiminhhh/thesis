print('1. Evaluate Unintegrated - Setup')
import os
N_CPU = 32
os.environ.setdefault("OMP_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
os.environ.setdefault("MKL_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
os.environ.setdefault("NUMEXPR_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
import scanpy as sc
from utils import score, align, eval_cells, joint_pca
from integration_datasets import DATASETS

method = 'Unintegrated'
for g, d in DATASETS.items():
    print(f'2. Evaluate {method} - {g}')
    raw = sc.read_h5ad(f"{d['concat']}.h5ad")
    adata = sc.read_h5ad(f"{d['concat']}.h5ad")
    adata, raw = align(adata, raw, eval_cells(d['out']))
    # Unintegrated embedding: PCA of the data scaled across both species
    joint_pca(adata)
    scib_scores = score(adata_int = adata,
                        method = method,
                        output_type = 'embed',
                        raw = raw,
                        embed = 'X_pca')

    print(f'3. Evaluate {method} - Save score - {g}')
    scib_scores.to_csv(f"{d['scores']}05.00_unintegrated.csv")
    print(f'Evaluate {method} - {g} - Done!')
