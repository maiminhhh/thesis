# Cells per SingleR (cross-species reference) label in one dataset, from the 02 h5ad; pruned calls -> 'Unassigned';
# cell types of the dataset's label set (same-species Seurat labels) that SingleR (cross) never assigns -> 0.
#   python -I analysis/singler_cross_counts.py <tissue> <species>  -> figures/<tissue>/annotation/singler_cross_counts_<species>.csv
import sys, os
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import anndata as ad
import pandas as pd

tissue, species = sys.argv[1:3]
obs = ad.read_h5ad(f'{ROOT}/{tissue}/data/02_annotated_data/02_{species}_annotated.h5ad', backed = 'r').obs
lab = obs['celltype_singleR_cross'].astype(object).where(obs['celltype_singleR_cross'].notna(), 'Unassigned')
n = lab.value_counts().rename_axis('label').reset_index(name = 'n_cells')
# reference types never assigned (label set: the same-species Seurat labels) shown with 0 cells
never = sorted(set(obs['celltype_seurat'].dropna().astype(str)) - set(n['label']))
n = pd.concat([n, pd.DataFrame({'label': never, 'n_cells': 0})], ignore_index = True)
n['pct'] = 100 * n['n_cells'] / n['n_cells'].sum()
out = f'{ROOT}/figures/{tissue}/annotation/singler_cross_counts_{species}.csv'
n.to_csv(out, index = False)
print(n.to_string(index = False)); print('total', n['n_cells'].sum(), '->', out)
