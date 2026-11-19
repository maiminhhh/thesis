print('1. scIB Setup')
print('Packages import')
import scanpy as sc
import pandas as pd
import scib

print('Methods info')
paths = {
    'unintegrated': 'data/03_concatenated_data.h5ad',
    'fastmnn': 'data/04_integrated_data/04.01_FastMNN.h5ad',
    'samap': 'data/04_integrated_data/04.02_SAMap.h5ad',
    'cca': 'data/04_integrated_data/04.03_CCA.h5ad',
    'harmony': 'data/04_integrated_data/04.04_Harmony.h5ad',
    'scanorama': 'data/04_integrated_data/04.05_Scanorama.h5ad',
    'rpca': 'data/04_integrated_data/04.06_RPCA.h5ad',
    'scvi': 'data/04_integrated_data/04.07_scVI.h5ad',
    'scanvi': 'data/04_integrated_data/04.08_scANVI.h5ad',
    'bbknn': 'data/04_integrated_data/04.09_bbknn.h5ad',
    'liger': 'data/04_integrated_data/04.10_LIGER.h5ad',
    'liger_uinmf': 'data/04_integrated_data/04.11_LIGER_UINMF.h5ad',
    'fastmnn_uinmf': 'data/04_integrated_data/04.12_UINMF_fastMNN.h5ad',
    'cca_uinmf': 'data/04_integrated_data/04.13_UINMF_CCA.h5ad',
    'harmony_uinmf': 'data/04_integrated_data/04.14_UINMF_Harmony.h5ad',
    'scanorama_uinmf': 'data/04_integrated_data/04.15_UINMF_Scanorama.h5ad',
    'rpca_uinmf': 'data/04_integrated_data/04.16_UINMF_RPCA.h5ad',
    'bbknn_uinmf': 'data/04_integrated_data/04.17_UINMF_BBKNN.h5ad'
}
embed_keys = {
    'unintegrated': 'X_pca',
    'fastmnn': 'X_mnn',
    'samap': 'wPCA',
    'cca': 'X_pca',
    'harmony': 'X_harmony',
    'scanorama': 'X_x_scanorama',
    'rpca': 'X_pca',
    'scvi': 'X_x_scvi',
    'scanvi': 'X_x_scanvi',
    'bbknn': 'X_pcabbknn',
    'liger': 'X_inmfnorm',
    'liger_uinmf': 'X_inmf',
    'fsatmnn_uinmf': 'X_fastmnn_uinmf',
    'cca_uinmf': 'X_uinmf_cca',
    'harmony_uinmf': 'X_harmony',
    'scanorama_uinmf': 'X_x_uinmf_scanorama',
    'rpca_uinmf': 'X_uinmf_rpca',
    'bbknn_uinmf': 'X_x_pca_bbknn_uinmf'
}