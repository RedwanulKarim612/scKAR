import pandas as pd
import os
import sys


def create_fasta(df, file_name):
    cnt = 0
    with open(file_name, 'w') as f:
        for index, row in df.iterrows():
            f.write('>kmer' + str(cnt) + '\n')
            f.write(index + '\n')
            cnt+=1


path = sys.argv[1]
log2fc = float(sys.argv[2])
padj = float(sys.argv[3])
base_mean = float(sys.argv[4])
bipartition_name = sys.argv[5]

deseq_folder_path = path + '/deseq_results/'
final_results_folder_path = path + '/final_results/' + bipartition_name + '/'

files = os.listdir(deseq_folder_path)
files = [f for f in files if f.endswith('.csv') and f.startswith(bipartition_name)]
merged_df = pd.DataFrame()

print('merging deseq results')

merged_df = pd.read_csv(deseq_folder_path + files[0], sep='\t', index_col=0)
filtered_df = merged_df.copy()

filtered_df = filtered_df[filtered_df['padj'] <= padj]
filtered_df = filtered_df[abs(filtered_df['log2FoldChange']) >= log2fc ]
filtered_df = filtered_df[filtered_df['baseMean'] >= base_mean]

# find Nan values
filtered_df.fillna(0, inplace=True)
filtered_df.index.name = 'kmer'
filtered_df.to_csv(final_results_folder_path + '/filtered_kmers.tsv', sep='\t')
print('created filtered_kmers.tsv')

print('creating fasta files')
A_df = filtered_df[filtered_df['log2FoldChange'] > 0]
B_df = filtered_df[filtered_df['log2FoldChange'] < 0]

create_fasta(A_df, final_results_folder_path + '/A_kmers.fasta')
create_fasta(B_df, final_results_folder_path + '/B_kmers.fasta')

A_df.to_csv(final_results_folder_path + '/A_kmers.tsv', sep='\t')
B_df.to_csv(final_results_folder_path + '/B_kmers.tsv', sep='\t')
print('created A_kmers.fasta and B_kmers.fasta')