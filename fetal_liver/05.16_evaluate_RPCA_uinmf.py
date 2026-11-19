print('1. Evaluate RPCA UINMF - Setup')
import os
N_CPU = 32
os.environ.setdefault("OMP_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
os.environ.setdefault("MKL_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
os.environ.setdefault("OPENBLAS_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
os.environ.setdefault("NUMEXPR_NUM_THREADS", str(N_CPU))   # run_05.sh caps threads
import scanpy as sc
from utils import score, align, eval_cells
from integration_datasets import DATASETS

method = 'RPCA UINMF'
for g, d in DATASETS.items():
    print(f'2. Evaluate {method} - {g}')
    raw = sc.read_h5ad(f"{d['concat']}.h5ad")
    adata = sc.read_h5ad(f"{d['out']}04.16_RPCA_UINMF.h5ad")
    adata, raw = align(adata, raw, eval_cells(d['out']))
    scib_scores = score(adata_int = adata,
                        method = method,
                        output_type = 'embed',
                        raw = raw,
                        embed = 'X_rpca_uinmf')

    print(f'3. Evaluate {method} - Save score - {g}')
    scib_scores.to_csv(f"{d['scores']}05.16_RPCA UINMF.csv")
    print(f'Evaluate {method} - {g} - Done!')
