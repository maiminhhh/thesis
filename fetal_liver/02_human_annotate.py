import tacco as tc
import scanpy as sc
import pandas as pd

ref = sc.read_h5ad('data/02_annotation_refs/human_ref.h5ad')
data = sc.read_h5ad('data/02_annotated_data/02_human_annotated.h5ad')
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

probs = data.obsm["celltype_tacco"]

valid = probs.notna().any(axis=1)

labels = pd.Series(index=probs.index, dtype="object")

labels.loc[valid] = probs.loc[valid].idxmax(axis=1)

data.obs["celltype_tacco"] = (
    labels.reindex(data.obs_names).astype("category")
)

data.write('data/02_annotated_data/02_human_annotated.h5ad')