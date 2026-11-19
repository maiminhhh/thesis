print('1. Evaluate LIGER UINMF - Setup')
import os
N_CPU = int(os.environ.get('N_THREADS', 32))
os.environ["OMP_NUM_THREADS"] = str(N_CPU)
os.environ["MKL_NUM_THREADS"] = str(N_CPU)
os.environ["OPENBLAS_NUM_THREADS"] = str(N_CPU)
os.environ["NUMEXPR_NUM_THREADS"] = str(N_CPU)
import scanpy as sc
from utils import score, align, eval_cells
from integration_datasets import DATASETS

method = 'LIGER UINMF'
for d, ds in DATASETS.items():
    print(f'2. Evaluate {method} - {d}')
    raw = sc.read_h5ad(f"{ds['concat']}.h5ad")
    adata = sc.read_h5ad(f"{ds['out']}04.11_LIGER_UINMF.h5ad")
    adata, raw = align(adata, raw, eval_cells(ds["out"]))
    scib_scores = score(adata_int = adata,
                        method = method,
                        output_type = 'embed',
                        raw = raw,
                        embed = 'X_inmf')

    print(f'3. Evaluate {method} - Save score - {d}')
    scib_scores.to_csv(f"{ds['scores']}05.11_LIGER UINMF.csv")
    print(f'Evaluate {method} - {d} - Done!')
