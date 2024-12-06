import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import sys

matrix_path     = sys.argv[1]
clustering_algo = sys.argv[2]
min_genes       = int(sys.argv[3])
min_cells       = int(sys.argv[4])
n_neighbors     = int(sys.argv[5])
n_pcs           = int(sys.argv[6])

adata = pd.read_csv(matrix_path, index_col=0)

x = adata.isna().sum()

# convert to scanpy object
adata = sc.AnnData(adata)
adata = adata.T
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_genes(adata, min_cells=2)
adata.obs['n_counts'] = adata.X.sum(axis=1)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.regress_out(adata, ['n_counts'], n_jobs=20)
sc.pp.scale(adata, max_value=10)
adata.obs['n_counts'] = adata.X.sum(axis=1)

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.1)
sc.pl.umap(adata, color='leiden')
# sc.tl.tsne(adata)
sc.tl.tsne(adata)
sc.pl.tsne(adata, color='leiden')
sc.tl.dendrogram(adata, groupby = 'leiden')
sc.pl.dendrogram(adata, groupby = 'leiden')

linkage_matrix = adata.uns['dendrogram_leiden']

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

print(child_set)
print(only_left_right_child_set)

# finding internal edges
universal_set = set(range(0, leaf_count))
bipartite_sets = []
for key in only_left_right_child_set: 
    for child in only_left_right_child_set[key]:
        if child >= leaf_count:
            # internal edge found
            bipartite_sets.append((child_set[child], universal_set.difference(child_set[child])))

print(bipartite_sets)

adata.obs['leiden'].to_csv('p64_2_leiden.csv')