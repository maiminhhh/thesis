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

print('4. Run SAMap')
# datasets (datasets.py): all by default, or e.g. `python 04.02_integrate_SAMap.py genes8` to run one dataset (datasets can run in parallel)
from datasets import selected
genes = selected()
for g in genes:
    human_adata = sc.read_h5ad(f'data/{g}/02_human_annotated.h5ad')
    mouse_subset_adata = sc.read_h5ad(f'data/{g}/02_mouse_annotated.h5ad')
    mouse_subset_adata.var_names = mouse_subset_adata.var['names']
    human_adata.var_names = human_adata.var['names']
    human = SAM(human_adata)
    mouse_subset = SAM(mouse_subset_adata)

    for a in (human.adata, mouse_subset.adata):
        for key in ['data', 'scale.data', 'X_disp', 'X_counts']:
            if key in a.layers:
                del a.layers[key]
                
    # prepend species prefix to gene names for SAMap
    prepend_var_prefix(human, 'hs')
    prepend_var_prefix(mouse_subset, 'mm')
    # preprocess data
    human.preprocess_data()
    human.run()
    mouse_subset.preprocess_data()
    mouse_subset.run()

    sams = {'hs': human, 'mm': mouse_subset}
    sm = SAMAP(sams,
            f_maps = 'data/04_integration_resources/maps/', backend = os.environ.get('SAMAP_BACKEND', 'cuda'), save_processed = True)   # SAMAP_BACKEND=cpu if the GPU runs out of memory
    sm.run()
    samap_adata = sm.samap.adata
    #samap_adata.obsm['wPCA'] = sm.samap.adata.obsm['X_umap'].copy()

    print('5. SAMap Save data')
    # 1. unified cell-type label: take mouse label for mm cells, human label for hs cells
    #    (pick whichever cell-type column you trust — using majority_vote here)
    hs_lab = sm.sams['hs'].adata.obs['celltype']   # indexed by hs barcodes
    mm_lab = sm.sams['mm'].adata.obs['celltype']   # indexed by mm barcodes

    sp = samap_adata.obs['species'].astype(str).values
    names = samap_adata.obs_names

    out = np.empty(samap_adata.n_obs, dtype=object)
    for i, (name, s) in enumerate(zip(names, sp)):
        src = mm_lab if s == 'mm' else hs_lab
        # strip any suffix SAMAP may have added to disambiguate duplicates
        key = name
        if key not in src.index and '-' in key:
            key = key.rsplit('-', 1)[0]          # try without trailing -N
        out[i] = src.get(key, np.nan)

    samap_adata.obs['celltype'] = out

    # 2. make species a clean string col (it already exists, just ensure dtype)
    samap_adata.obs['orig.ident'] = samap_adata.obs['species'].astype(str)

    # 3. cast ALL obs columns to str so the R converter doesn't choke on category/NaN
    for c in samap_adata.obs.columns:
        samap_adata.obs[c] = samap_adata.obs[c].astype(str)

    # 4. keep real gene names (mm_/hs_) as the var index
    samap_adata.var_names = samap_adata.var_names.astype(str)

    # 5. sanity check before writing
    print(samap_adata.obs[['orig.ident','celltype']].head())
    print(samap_adata.obs['orig.ident'].value_counts())
    print(samap_adata.var_names[:5].tolist(), samap_adata.var_names[-5:].tolist())

    print('5. SAMap Save data')
    samap_adata.write(f'data/{g}/04.02_SAMap.h5ad')


    print('SAMap Done!')