# md5 of the cell and gene names of every 03 concatenated dataset (section 5 of ../EBCC_PIPELINE.md)
# usage: python fingerprint_03.py > out.txt
import hashlib, warnings
import anndata as ad
warnings.filterwarnings('ignore')
files = {'full_mouse': 'data/03_concatenated_data/03_concatenated_data_full.h5ad',
         'dissected_mouse': 'data/03_concatenated_data/03_concatenated_data_subset.h5ad',
         **{f'genes{g}': f'data/genes{g}/03_concatenated_data.h5ad' for g in (2, 4, 6, 8)}}
md5 = lambda x: hashlib.md5('\n'.join(map(str, x)).encode()).hexdigest()
for name, f in files.items():
    a = ad.read_h5ad(f, backed='r')
    print(name, a.n_obs, a.n_vars, md5(a.obs_names), md5(a.var_names))
