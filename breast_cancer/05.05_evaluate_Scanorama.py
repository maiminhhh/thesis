print('1. Evaluate Scanorama - Setup')
import os
N_CPU = int(os.environ.get('N_CPU', 32))   # threads; the runner sets this
os.environ["OMP_NUM_THREADS"] = str(N_CPU)
os.environ["MKL_NUM_THREADS"] = str(N_CPU)
os.environ["OPENBLAS_NUM_THREADS"] = str(N_CPU)
os.environ["NUMEXPR_NUM_THREADS"] = str(N_CPU)
import scanpy as sc
from utils import score, align

method = 'Scanorama'
# datasets (datasets.py): all by default, or e.g. `python 05.05_evaluate_Scanorama.py genes8` to score one dataset
from datasets import selected, scores_dir
genes = selected()
for g in genes:
    print(f'2. Evaluate {method} - {g}')
    raw = sc.read_h5ad(f'data/{g}/03_concatenated_data.h5ad')
    adata = sc.read_h5ad(f'data/{g}/04.05_Scanorama.h5ad')
    adata, raw = align(adata, raw, g)
    scib_scores = score(adata_int = adata,
                        method = method,
                        output_type = 'embed',
                        raw = raw,
                        embed = 'X_scanorama')

    print(f'3. Evaluate {method} - Save score - {g}')
    scib_scores.to_csv(scores_dir(g) + '05.05_Scanorama.csv')
    print(f'Evaluate {method} - {g} - Done!')
