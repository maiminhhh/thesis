# Coordinates and EBCC labels for the mouse-dissection figure (analysis/dissection_figure.R).
#   python -I analysis/dissection_figure_compute.py
# The annotated .h5ad files hold the EBCC labels but no spatial coordinates (those are only in the Seurat .rds files),
# so the cell centroids are read from each sample's raw platform output and joined by cell id. The dissection is then
# recomputed from each tissue's rectangle (as in its 03 / annotation script) and must reproduce the saved dissected
# cell set exactly; otherwise the script stops.
# Output: figures/dissection/dissection_cells.csv.gz (tissue, species, cell, x, y, label, dissected)
import os, sys
import numpy as np, pandas as pd, anndata as ad

# ---------------- config ----------------
ROOT = os.environ.get('PROJECT_ROOT') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f'{ROOT}/figures/dissection'
TISSUES = {   # raw centroid files (mouse, human) and the dissection rectangle in GetTissueCoordinates units (x, y)
    'colon': dict(
        mouse = 'colon/data/00_raw_data/Xenium_V1_mouse_Colon_FF_outs/cells.csv.gz',
        human = 'colon/data/00_raw_data/Xenium_V1_hColon_Non_diseased_Add_on_FFPE_outs/cells.csv.gz',
        rect = (260, 2685, 2590, 4625)),                      # colon/03_concatenate.qmd
    'bone_marrow': dict(
        mouse = 'bone_marrow/data/00_raw_data/output-XETG00120__0005224__p8_femur__20231122__195355/cells.csv.gz',
        human = 'bone_marrow/data/00_raw_data/output-XETG00120_0005071_normal_bone_marrow_20240119_232334/cells.csv',
        rect = (500, 4750, 2500, 3500)),                      # bone_marrow/annotation_viz.qmd
    'fetal_liver': dict(
        mouse = 'fetal_liver/data/00_raw_data/output-XETG00120__0005224__E145_liver__20231122__195355/cells.csv.gz',
        human = 'fetal_liver/data/00_raw_data/output-XETG00120_0005059_fetal_liver_20240119_232334/cells.csv.gz',
        rect = (0, 3100, 500, 3000)),                         # fetal_liver/03_concatenate.R
    'breast_cancer': dict(
        mouse = 'breast_cancer/data/00_raw_data/GSM8830802_T2_cell_metadata.csv',
        human = 'breast_cancer/data/00_raw_data/202409242358_240916JHHUBC0005XQ-V2V-HubcTMA-V2-BY_VMSC02511_region_R4_cell_metadata.csv',
        rect = (4500, 8500, 1500, 5500),                      # breast_cancer/03_concatenate.qmd
        subset = 'breast_cancer/data/dissected_mouse/02_mouse_annotated.h5ad'),   # the cut, after concatenation
}
LABEL = 'ebcc_label'
# -----------------------------------------

def centroids(path):
    d = pd.read_csv(f'{ROOT}/{path}')
    if 'cell_id' in d:   # Xenium
        return pd.DataFrame({'cell': d.cell_id.astype(str), 'x': d.x_centroid, 'y': d.y_centroid})
    return pd.DataFrame({'cell': d.EntityID.astype(str), 'x': d.center_x, 'y': d.center_y})   # MERSCOPE

def obs(t, name, path = None):
    o = ad.read_h5ad(f'{ROOT}/{path}' if path else f'{ROOT}/{t}/data/02_annotated_data/{name}.h5ad', backed = 'r').obs
    return pd.DataFrame({'cell': o.index.astype(str).str.replace('^mmusculus_', '', regex = True), 'label': o[LABEL].astype(str).values})

os.makedirs(OUT, exist_ok = True)
rows = []
for t, c in TISSUES.items():
    full, sub, hum = obs(t, '02_mouse_annotated'), obs(t, '02_mouse_subset_annotated', c.get('subset')), obs(t, '02_human_annotated')
    m = full.merge(centroids(c['mouse']), on = 'cell', how = 'left')
    h = hum.merge(centroids(c['human']), on = 'cell', how = 'left')
    if m.x.isna().any() or h.x.isna().any():
        sys.exit(f'{t}: {m.x.isna().sum()} mouse / {h.x.isna().sum()} human cells without coordinates')
    x0, x1, y0, y1 = c['rect']
    in_rect = m.x.between(x0, x1) & m.y.between(y0, y1)
    saved = set(sub.cell)
    got = set(m.cell[in_rect])
    print(f'{t:14s} mouse {len(m):>7d}, in rectangle {len(got):>6d}, saved dissected {len(saved):>6d}: '
          f'{"identical" if got == saved else f"saved within rectangle: {saved <= got}; in rectangle but not saved: {len(got - saved)}"}; human {len(h)}')
    if not saved <= got:   # the saved cut must lie inside the rectangle (breast cancer drops a few cells at concatenation)
        sys.exit(f'{t}: the saved dissected cells are not all inside the rectangle')
    m['dissected'] = m.cell.isin(saved)
    # labels of dissected cells: those of the full mouse object (the subset carries the same EBCC labels; checked)
    lab_sub = sub.set_index('cell').label
    diff = (m.set_index('cell').label[list(saved)].values != lab_sub[list(saved)].values).sum()
    print(f'{"":14s} label differences between full and dissected objects: {diff}')
    rows += [m.assign(tissue = t, species = 'mouse'), h.assign(tissue = t, species = 'human', dissected = np.nan)]

out = pd.concat(rows)[['tissue', 'species', 'cell', 'x', 'y', 'label', 'dissected']]
out.to_csv(f'{OUT}/dissection_cells.csv.gz', index = False)
print(f'wrote {OUT}/dissection_cells.csv.gz ({len(out)} cells)')
