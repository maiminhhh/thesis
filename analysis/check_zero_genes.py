# Panel genes with zero detection: for each tissue and species, every feature in the raw cell-by-gene matrix
# (all cells the instrument segmented, before any QC), total counts and the number of cells detecting it.
# Xenium: cell_feature_matrix.h5 (feature_type separates gene targets from control probes/codewords);
# breast cancer (Vizgen MERSCOPE): cell_by_gene.csv (Blank-* columns are controls). Reads raw files only.
# Usage: python analysis/check_zero_genes.py  -> figures/cross_tissue/raw_gene_detection.csv (all features), printed summary
import h5py, numpy as np, pandas as pd, scipy.sparse as sp

import os
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XENIUM = {
    ('colon', 'human'): 'colon/data/00_raw_data/Xenium_V1_hColon_Non_diseased_Add_on_FFPE_outs',
    ('colon', 'mouse'): 'colon/data/00_raw_data/Xenium_V1_mouse_Colon_FF_outs',
    ('bone_marrow', 'human'): 'bone_marrow/data/00_raw_data/output-XETG00120_0005071_normal_bone_marrow_20240119_232334',
    ('bone_marrow', 'mouse'): 'bone_marrow/data/00_raw_data/output-XETG00120__0005224__p8_femur__20231122__195355',
    ('fetal_liver', 'human'): 'fetal_liver/data/00_raw_data/output-XETG00120_0005059_fetal_liver_20240119_232334',
    ('fetal_liver', 'mouse'): 'fetal_liver/data/00_raw_data/output-XETG00120__0005224__E145_liver__20231122__195355',
}
MERSCOPE = {
    ('breast_cancer', 'human'): 'breast_cancer/data/00_raw_data/202409242358_240916JHHUBC0005XQ-V2V-HubcTMA-V2-BY_VMSC02511_region_R4_cell_by_gene.csv',
    ('breast_cancer', 'mouse'): 'breast_cancer/data/00_raw_data/GSM8830802_T2_cell_by_gene.csv',
}

rows = []
for (t, s), d in XENIUM.items():
    with h5py.File(f'{ROOT}/{d}/cell_feature_matrix.h5') as f:
        m = f['matrix']
        n_feat, n_cell = m['shape'][:]
        X = sp.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape = (n_feat, n_cell)).tocsr()   # features x cells
        names = [x.decode() for x in m['features/name'][:]]
        ftype = [x.decode() for x in m['features/feature_type'][:]]
    rows.append(pd.DataFrame({'tissue': t, 'species': s, 'feature': names, 'feature_type': ftype, 'n_cells_total': n_cell,
                              'total_counts': np.asarray(X.sum(1)).ravel(), 'n_cells_detected': np.diff(X.indptr)}))
for (t, s), f in MERSCOPE.items():
    df = pd.read_csv(f'{ROOT}/{f}', index_col = 0)
    v = df.to_numpy()
    rows.append(pd.DataFrame({'tissue': t, 'species': s, 'feature': df.columns,
                              'feature_type': np.where(df.columns.str.startswith('Blank'), 'Blank Codeword', 'Gene Expression'),
                              'n_cells_total': v.shape[0], 'total_counts': v.sum(0), 'n_cells_detected': (v > 0).sum(0)}))
d = pd.concat(rows, ignore_index = True)
d['pct_cells_detected'] = 100 * d.n_cells_detected / d.n_cells_total
d.to_csv(f'{ROOT}/figures/cross_tissue/raw_gene_detection.csv', index = False)

g = d[d.feature_type == 'Gene Expression']
print('feature types:', d.groupby('feature_type').size().to_dict(), '\n')
summ = g.groupby(['tissue', 'species']).agg(n_cells = ('n_cells_total', 'first'), n_genes = ('feature', 'size'),
                                            zero = ('total_counts', lambda x: int((x == 0).sum())),
                                            under_0_1pct = ('pct_cells_detected', lambda x: int((x < 0.1).sum())),
                                            under_1pct = ('pct_cells_detected', lambda x: int((x < 1).sum())),
                                            median_pct_cells = ('pct_cells_detected', 'median'))
print(summ.round(2).to_string(), '\n')
print('genes with zero counts in all cells:')
print(g[g.total_counts == 0][['tissue', 'species', 'feature', 'n_cells_total']].to_string(index = False) or 'none')
print('\nleast-detected genes per dataset:')
print(g.sort_values('pct_cells_detected').groupby(['tissue', 'species']).head(3)
       [['tissue', 'species', 'feature', 'total_counts', 'n_cells_detected', 'pct_cells_detected']].round(4).to_string(index = False))
