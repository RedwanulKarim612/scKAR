import pandas as pd
import sys

tpm_sum_path = sys.argv[1]
clustering_results_path = sys.argv[2]
output_path = sys.argv[3]
tpm_sum = pd.read_csv(tpm_sum_path, index_col = 0, header=None, sep='\t')
tpm_sum.columns = ['tpm_sum']

clustering_results = pd.read_csv(clustering_results_path, index_col = 0)

tpm_cluster_df = clustering_results.merge(tpm_sum, left_index = True, right_index = True)
tpm_cluster_df.sort_index(inplace = True)
tpm_cluster_df.to_csv(sys.argv[3], index = True, columns = ['tpm_sum', 'cluster'], header = None, sep = '\t')