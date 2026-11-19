print('1. SAMap Setup')
import pandas as pd
import scanpy as sc
import numpy as np
import anndata as ad
import samap
from samap import SAMAP
from samap import get_mapping_scores, sankey_plot
from samap.utils import prepend_var_prefix
from samap.sam import SAM
import re
import gzip
from collections import defaultdict
import pickle
import os
from datasets import selected

print('4. Run SAMap')
for name, d in selected():
    print(f'SAMap - {name}')
    human_adata = sc.read_h5ad(d['human'])
    human_adata.var_names = human_adata.var['names']
    mouse_subset_adata = sc.read_h5ad(d['mouse'])
    mouse_subset_adata.var_names = mouse_subset_adata.var['names']
    human = SAM(human_adata)
    mouse_subset = SAM(mouse_subset_adata)

    for a in (human.adata, mouse_subset.adata):
        for key in ['data', 'scale.data', 'X_disp', 'X_counts']:
            if key in a.layers:
                del a.layers[key]

    prepend_var_prefix(human, 'hs')
    prepend_var_prefix(mouse_subset, 'mm')

    human.preprocess_data()
    human.run()
    mouse_subset.preprocess_data()
    mouse_subset.run()

    sams = {'hs': human, 'mm': mouse_subset}
    sm = SAMAP(sams,
            f_maps = 'data/04_integration_resources/maps/', backend = 'cuda', save_processed = True)
    sm.run()
    samap_adata = sm.samap.adata

    print('5. SAMap Save data')
    hs_lab = sm.sams['hs'].adata.obs['ebcc_label']   # indexed by hs barcodes
    mm_lab = sm.sams['mm'].adata.obs['ebcc_label']   # indexed by mm barcodes

    sp = samap_adata.obs['species'].astype(str).values
    names = samap_adata.obs_names

    out = np.empty(samap_adata.n_obs, dtype=object)
    for i, (cell, s) in enumerate(zip(names, sp)):
        src = mm_lab if s == 'mm' else hs_lab
        key = cell
        if key not in src.index and '-' in key:
            key = key.rsplit('-', 1)[0]          # try without trailing -N
        out[i] = src.get(key, np.nan)

    samap_adata.obs['ebcc_label'] = out
    samap_adata.obs['orig.ident'] = samap_adata.obs['species'].astype(str)

    for c in samap_adata.obs.columns:
        samap_adata.obs[c] = samap_adata.obs[c].astype(str)

    samap_adata.var_names = samap_adata.var_names.astype(str)

    print(samap_adata.obs[['orig.ident','ebcc_label']].head())
    print(samap_adata.obs['orig.ident'].value_counts())
    print(samap_adata.var_names[:5].tolist(), samap_adata.var_names[-5:].tolist())

    samap_adata.write(f"{d['out']}04.02_SAMap.h5ad")
    print(f'SAMap Done! - {name}')
