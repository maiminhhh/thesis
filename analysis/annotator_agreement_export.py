# Cross-tab of two annotator columns of one dataset, from the 02 h5ad (as annotation_agreement() in analysis/visualisations.R,
# which needs the Seurat .rds): rows = <row_col> labels, columns = <col_col>; missing (SingleR pruned) -> 'Unassigned';
# frac = share of the row label's cells.
#   python -I analysis/annotator_agreement_export.py <tissue> <species> <row_col> <col_col> <out.csv>
# Output columns: seurat (row label), singleR (column label), n, frac, species  (names as plot_annotation_agreement expects)
import sys, os
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import anndata as ad
import pandas as pd

tissue, species, row_col, col_col, out = sys.argv[1:6]
obs = ad.read_h5ad(f'{ROOT}/{tissue}/data/02_annotated_data/02_{species}_annotated.h5ad', backed = 'r').obs
lab = lambda c: obs[c].astype(object).where(obs[c].notna(), 'Unassigned').astype(str)
d = pd.DataFrame({'seurat': lab(row_col), 'singleR': lab(col_col)}).value_counts().reset_index(name = 'n')
d['frac'] = d['n'] / d.groupby('seurat')['n'].transform('sum')
d['species'] = species
d.sort_values(['seurat', 'singleR']).to_csv(out, index = False)
print(f'{tissue} {species}: {row_col} x {col_col}, {d.n.sum()} cells -> {out}')
