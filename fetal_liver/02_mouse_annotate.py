import tacco as tc
import scanpy as sc
import pandas as pd

ref = sc.read_h5ad('data/02_annotation_refs/02_mouse_ref.h5ad')
data = sc.read_h5ad('data/02_annotated_data/02_mouse_annotated.h5ad')
ref = tc.pp.construct_reference_profiles(
    ref,
    annotation_key='celltype',
    inplace=False
)

ref.varm['celltype'] = ref.varm['celltype'].copy()

tc.tl.annotate(
    data,
    ref,
    annotation_key='celltype',
    result_key='celltype_tacco'
)

data.varm['celltype_tacco'].columns = (
    data.varm['celltype_tacco'].columns
    .str.replace('/', '_', regex=False)
)

data.obsm['celltype_tacco'].columns = (
    data.obsm['celltype_tacco'].columns
    .str.replace('/', '_', regex=False)
)
probs = data.obsm["celltype_tacco"]

# Find rows with at least one valid probability
valid = probs.notna().any(axis=1)

# Initialise labels as missing
labels = pd.Series(index=probs.index, dtype="object")

# Assign the most probable cell type to valid rows
labels.loc[valid] = probs.loc[valid].idxmax(axis=1)

# Align labels to AnnData observation order
data.obs["celltype_tacco"] = (
    labels.reindex(data.obs_names).astype("category")
)

data.write('data/02_annotated_data/02_mouse_annotated.h5ad')