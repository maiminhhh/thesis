print('1. Evaluate Scanorama scVI-latent - Setup')
import os
N_CPU = 32
os.environ.setdefault("OMP_NUM_THREADS", str(N_CPU))
os.environ.setdefault("MKL_NUM_THREADS", str(N_CPU))
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(N_CPU))
os.environ.setdefault("NUMEXPR_NUM_THREADS", str(N_CPU))
import scanpy as sc
from utils import score, align, eval_cells
from datasets import selected

method = 'Scanorama scVI-latent'
for name, d in selected():
    print(f'2. Evaluate {method} - {name}')
    raw = sc.read_h5ad(d['concat'])
    adata = sc.read_h5ad(f"{d['out']}04.22_Scanorama_scVI.h5ad")
    adata, raw = align(adata, raw, eval_cells(d['out']))
    scib_scores = score(adata_int = adata,
                        method = method,
                        output_type = 'embed',
                        raw = raw,
                        embed = 'X_scanorama_scvi')

    print(f'3. Evaluate {method} - Save score - {name}')
    scib_scores.to_csv(f"{d['scores']}05.22_Scanorama scVI-latent.csv")
    print(f'Evaluate {method} - {name} - Done!')
