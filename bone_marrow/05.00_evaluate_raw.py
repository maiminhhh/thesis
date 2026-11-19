print('1. Evaluate Unintegrated - Setup')
import os
N_CPU = 32
os.environ.setdefault("OMP_NUM_THREADS", str(N_CPU))
os.environ.setdefault("MKL_NUM_THREADS", str(N_CPU))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(N_CPU))
os.environ.setdefault("NUMEXPR_NUM_THREADS", str(N_CPU))
import scanpy as sc
from utils import score, align, eval_cells, joint_pca
from datasets import selected

method = 'Unintegrated'
for name, d in selected():
    print(f'2. Evaluate {method} - {name}')
    raw = sc.read_h5ad(d['concat'])
    adata = sc.read_h5ad(d['concat'])
    adata, raw = align(adata, raw, eval_cells(d['out']))
    joint_pca(adata)   # same jointly scaled PCA as the PCR baseline
    scib_scores = score(adata_int = adata,
                        method = method,
                        output_type = 'embed',
                        raw = raw,
                        embed = 'X_pca')

    print(f'3. Evaluate {method} - Save score - {name}')
    scib_scores.to_csv(f"{d['scores']}05.00_unintegrated.csv")
    print(f'Evaluate {method} - {name} - Done!')
