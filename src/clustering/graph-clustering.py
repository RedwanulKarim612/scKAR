import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

matrix_path           = sys.argv[1]
clustering_algo       = sys.argv[2]        # possible options: leiden, louvain
min_genes             = int(sys.argv[3])   # by default 10
min_cells             = int(sys.argv[4])   # by default 2
n_neighbors           = int(sys.argv[5])   # by default 10
n_pcs                 = int(sys.argv[6])   # by default 40
clustering_resolution = float(sys.argv[7]) # by default 0.4

# get directory name of matrix path
directory = matrix_path.split('/')
directory = '/'.join(directory[:-2])

if not os.path.exists(directory + '/clustering_results'):
    os.mkdir(directory + '/clustering_results')

clustering_results_dir = directory + '/clustering_results/'

sc.settings.figdir = clustering_results_dir + '/figures'

if clustering_algo not in ['leiden', 'louvain']:
    print('clustering_algo must be either leiden or louvain')
    sys.exit(1)

adata = pd.read_csv(matrix_path, index_col=0)

x = adata.isna().sum()

# convert to scanpy object
adata = sc.AnnData(adata)
adata = adata.T
sc.pp.filter_cells(adata, min_genes = min_genes)
sc.pp.filter_genes(adata, min_cells = min_cells)
adata.obs['n_counts'] = adata.X.sum(axis=1)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.regress_out(adata, ['n_counts'], n_jobs=20)
sc.pp.scale(adata, max_value=10)
adata.obs['n_counts'] = adata.X.sum(axis=1)

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors = n_neighbors, n_pcs = n_pcs)

if clustering_algo == 'leiden':
    sc.tl.leiden(adata, resolution = clustering_resolution)
else:
    sc.tl.louvain(adata, resolution = clustering_resolution)
sc.tl.umap(adata)
sc.pl.umap(adata, color=clustering_algo, save =  '-' + clustering_algo + '.png')

sc.tl.tsne(adata)
sc.pl.tsne(adata, color=clustering_algo)

sc.tl.dendrogram(adata, groupby = clustering_algo)
sc.pl.dendrogram(adata, groupby = clustering_algo, save = '-' + clustering_algo + '.png')

linkage_matrix = adata.uns['dendrogram_' + clustering_algo]

# creating child list
leaf_count = len(linkage_matrix['dendrogram_info']['leaves'])
linkage = linkage_matrix['linkage']
child_set = {}
only_left_right_child_set = {}
for i in range(leaf_count):
    child_set[i] = {i}
cur_node = leaf_count
for link in linkage:
    node1 = int(link[0])
    node2 = int(link[1])
    only_left_right_child_set[cur_node] = {node1, node2}
    child_set[cur_node] = child_set[node1].union(child_set[node2])
    cur_node += 1

# finding internal edges
universal_set = set(range(0, leaf_count))
bipartite_sets = []
for key in only_left_right_child_set: 
    for child in only_left_right_child_set[key]:
        if child >= leaf_count:
            # internal edge found
            bipartite_sets.append((child_set[child], universal_set.difference(child_set[child])))


cluster_df = adata.obs[clustering_algo]
cluster_df = cluster_df.rename('cluster')
# name index column as cell
cluster_df.index.name = 'cell'
cluster_df.to_csv(clustering_results_dir + '/cluster.csv')

if not os.path.exists(clustering_results_dir + '/bipartitions'):
    os.mkdir(clustering_results_dir + '/bipartitions')
for i in range(len(bipartite_sets)):
    temp = cluster_df.copy()
    temp = temp.reset_index()

    temp['cluster'] = temp['cluster'].replace(list(str(bipartite_sets[i][0])), 'A')
    temp['cluster'] = temp['cluster'].replace(list(str(bipartite_sets[i][1])), 'B')
    temp.columns = ['cell', 'condition']

    temp.sort_values(by=['cell'], inplace=True)
    temp.to_csv(clustering_results_dir + "/bipartitions/bipartition_" + str(i) + ".csv", columns=["cell", "condition"], index=False, header=False)
