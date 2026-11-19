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
from collections import defaultdict
import os
from integration_datasets import DATASETS

for name, d in DATASETS.items():
    print(f'===== {name} =====')
    out_path = d['out']
    print('2. SAMap Preprocessing')
    # per-species inputs from 03: Ensembl ids (what the BLAST maps use), cells matched to the concat
    human_adata = sc.read_h5ad(f"{d['inputs']}02_human_annotated.h5ad")
    mouse_adata = sc.read_h5ad(f"{d['inputs']}02_mouse_annotated.h5ad")
    assert human_adata.var_names[0].startswith('ENSG') and mouse_adata.var_names[0].startswith('ENSMUSG'), \
        'SAMap inputs must use Ensembl ids'

    human = SAM(human_adata)
    mouse = SAM(mouse_adata)

    for a in (human.adata, mouse.adata):
        for key in ['data', 'scale.data', 'X_disp', 'X_counts']:
            if key in a.layers:
                del a.layers[key]
                
    # prepend species prefix to gene names for SAMap
    prepend_var_prefix(human, 'hs')
    prepend_var_prefix(mouse, 'mm')
    # preprocess data
    human.preprocess_data()
    human.run()
    mouse.preprocess_data()
    mouse.run()

    print('3. SAMap BLAST')
    # def fasta_to_id(input_fa, output_fa, prefix, keep_genes=None):
    #     """
    #     Rewrites FASTA headers from protein/transcript IDs to prefixed Ensembl gene IDs.
    #     keep_genes: set of *unprefixed, unversioned* gene IDs to keep (e.g. {'ENSMUSG00000021135', ...})
    #     """
    #     seen_genes = set()
    #     written = 0
    #     skipped = 0
    #     with open(input_fa, 'rt') as fin, open(output_fa, 'w') as fout:
    #         write_seq = False
    #         for line in fin:
    #             if line.startswith('>'):
    #                 gene_match = re.search(r'gene:(\S+)', line)
    #                 if gene_match:
    #                     gid = gene_match.group(1).split('.')[0]   # strip version suffix
    #                     prefixed = f"{prefix}_{gid}"
    #                     # only keep genes in the data, skip duplicates (one protein per gene)
    #                     if (keep_genes is None or gid in keep_genes) and prefixed not in seen_genes:
    #                         fout.write(f'>{prefixed}\n')
    #                         seen_genes.add(prefixed)
    #                         write_seq = True
    #                         written += 1
    #                     else:
    #                         write_seq = False
    #                         skipped += 1
    #                 else:
    #                     write_seq = False
    #             elif write_seq:
    #                 fout.write(line)
    #     print(f"Written: {written} genes, Skipped: {skipped}")
    #     return seen_genes


    # # get your gene IDs without prefix or version
    # human_genes = {g.replace('hs_', '').split('.')[0] for g in human.adata.var_names}
    # mouse_genes = {g.replace('mm_', '').split('.')[0] for g in mouse.adata.var_names}

    # print(f"Human genes to match: {len(human_genes)}")
    # print(f"Mouse genes to match: {len(mouse_genes)}")


    # DATA = 'data/04_integration_resources'

    # human_fa_out = f'{DATA}/hs_ids.fa'
    # mouse_fa_out = f'{DATA}/mm_ids.fa'

    # # --- sanity check overlap BEFORE the full rewrite ---
    # def check_overlap(input_fa, keep_genes):
    #     fa_ids = set()
    #     with open(input_fa) as f:
    #         for line in f:
    #             if line.startswith('>'):
    #                 m = re.search(r'gene:(\S+)', line)
    #                 if m:
    #                     fa_ids.add(m.group(1).split('.')[0])
    #     print(f"  overlap: {len(fa_ids & keep_genes)} of {len(keep_genes)}")

    # print("Human:"); check_overlap(f'{DATA}/Homo_sapiens.GRCh38.pep.all.fa', human_genes)
    # print("Mouse:"); check_overlap(f'{DATA}/Mus_musculus.GRCm39.pep.all.fa', mouse_genes)

    # # --- if the overlap looks right (near your panel size), run the rewrite ---
    # fasta_to_id(
    #     input_fa   = f'{DATA}/Homo_sapiens.GRCh38.pep.all.fa',
    #     output_fa  = human_fa_out,
    #     prefix     = 'hs',
    #     keep_genes = human_genes,
    # )

    # fasta_to_id(
    #     input_fa   = f'{DATA}/Mus_musculus.GRCm39.pep.all.fa',
    #     output_fa  = mouse_fa_out,
    #     prefix     = 'mm',
    #     keep_genes = mouse_genes,
    # )

    # %%bash
    # set -e

    # DATA=data/04_integration_resources
    # HS_FA=$DATA/hs_ids.fa
    # MM_FA=$DATA/mm_ids.fa

    # # make output dirs
    # mkdir -p $DATA/maps/mmhs
    # mkdir -p $DATA/maps/hsmm

    # # make blast databases
    # makeblastdb -in $MM_FA -dbtype prot -out $DATA/maps/mm_db
    # makeblastdb -in $HS_FA -dbtype prot -out $DATA/maps/hs_db

    # # blast mm -> hs
    # blastp -query $MM_FA \
    #        -db $DATA/maps/hs_db \
    #        -outfmt 6 \
    #        -out $DATA/maps/mmhs/mm_to_hs.txt \
    #        -num_threads 8 \
    #        -evalue 1e-6

    # # blast hs -> mm
    # blastp -query $HS_FA \
    #        -db $DATA/maps/mm_db \
    #        -outfmt 6 \
    #        -out $DATA/maps/hsmm/hs_to_mm.txt \
    #        -num_threads 8 \
    #        -evalue 1e-6

    # echo "BLAST done!"

    print('4. Run SAMap')
    sams = {'hs': human, 'mm': mouse}
    sm = SAMAP(sams,
            f_maps = 'data/04_integration_resources/maps/', backend = os.environ.get('SAMAP_BACKEND', 'cuda'), save_processed = True)
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
    samap_adata.write(f'{out_path}04.02_SAMap.h5ad')

    print('SAMap Done!')