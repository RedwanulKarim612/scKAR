import pandas as pd
import sys
import os

input_dir = sys.argv[1]
clustering_results_dir = input_dir + '/clustering_results'

if not os.path.exists(clustering_results_dir):
    os.mkdir(clustering_results_dir)
if not os.path.exists(clustering_results_dir + '/bipartitions'):
    os.mkdir(clustering_results_dir + '/bipartitions')

cluster_df = pd.read_csv(input_dir + '/metadata/clustering.csv', index_col=None, header=None)
cluster_df.columns = ['cell', 'cluster']

bipartitions = pd.read_csv(input_dir + '/metadata/bipartitions.csv', index_col=None, sep='\t')

for index, row in bipartitions.iterrows():
    set1 = eval(row['set1'])
    set2 = eval(row['set2'])

    temp = cluster_df.copy()
    temp['cluster'] = temp['cluster'].apply(lambda x: 'A' if x in set1 else ('B' if x in set2 else x))
    temp.sort_values(by=['cell'], inplace=True)
    temp.to_csv(clustering_results_dir + '/bipartitions/bipartition_' + str(index) + '.csv',
        sep='\t',
        index=False,
        header=False
    )

# copy the provided clusters to the clustering_results directory
cluster_df.to_csv(clustering_results_dir + '/cluster.csv', index=False)