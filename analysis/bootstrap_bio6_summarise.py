# The analysis/bootstrap_scores.py summaries with bio = all six bio metrics (as the benchmark: analysis/visualisations.R,
# EBCC_PIPELINE.md), instead of the four in bootstrap_scores.BIO_METRICS. The raw scores already hold the
# isolated-label metrics, so nothing is rescored: bootstrap_bio6/ links to the bootstrap/ score folders (read only)
# and bootstrap_scores.summarise() runs unchanged with the longer metric list, writing bootstrap_bio6/summary/.
# Usage: python analysis/bootstrap_bio6_summarise.py
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bootstrap_scores as bs
bs.OUT = f'{bs.ROOT}/bootstrap_bio6'
bs.BIO_METRICS = ['NMI_cluster/label', 'ARI_cluster/label', 'ASW_label',
                  'isolated_label_F1', 'isolated_label_silhouette', 'cLISI']
bs.summarise()
