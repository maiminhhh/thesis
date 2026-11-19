print('1. Evaluate BBKNN - Setup')
import os
N_CPU = 32
os.environ.setdefault("OMP_NUM_THREADS", str(N_CPU))
os.environ.setdefault("MKL_NUM_THREADS", str(N_CPU))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(N_CPU))
os.environ.setdefault("NUMEXPR_NUM_THREADS", str(N_CPU))
import scanpy as sc
from utils import score, align, eval_cells
from datasets import selected

method = 'BBKNN'
for name, d in selected():
    print(f'2. Evaluate {method} - {name}')
    raw = sc.read_h5ad(d['concat'])
    adata = sc.read_h5ad(f"{d['out']}04.09_BBKNN.h5ad")
    adata, raw = align(adata, raw, eval_cells(d['out']))
    scib_scores = score(adata_int = adata,
                        method = method,
                        output_type = 'knn',
                        raw = raw)

    print(f'3. Evaluate {method} - Save score - {name}')
    scib_scores.to_csv(f"{d['scores']}05.09_BBKNN.csv")
    print(f'Evaluate {method} - {name} - Done!')
