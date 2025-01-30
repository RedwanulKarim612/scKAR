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
files = os.listdir(path)
files = [f for f in files if f.endswith('.csv')]
print(len(files))
merged_df = pd.DataFrame()

print('merging deseq results')
# for file in files:
#     tmp_df = pd.read_csv(path + file, sep='\t')
#     merged_df = pd.concat([merged_df, tmp_df], axis=0)
    # print(file, merged_df.shape)

merged_df = pd.read_csv(path + files[0], sep='\t', index_col=0)
# merged_df = merged_df.iloc[:, :6]
filtered_df = merged_df.copy()
filtered_df = filtered_df[filtered_df['padj'] <= 0.05]

# filtered_df = filtered_df[abs(filtered_df['log2FoldChange']) >= 1 ]
# filtered_df = filtered_df[filtered_df['baseMean'] >= 8.0]
# find Nan values
filtered_df.fillna(0, inplace=True)
filtered_df.index.name = 'kmer'
filtered_df.to_csv(path + 'filtered_kmers.tsv', sep='\t')
print('created filtered_kmers.tsv')

print('creating fasta files')
A_df = filtered_df[filtered_df['log2FoldChange'] > 0]
B_df = filtered_df[filtered_df['log2FoldChange'] < 0]

create_fasta(A_df, path + 'A_kmers.fasta')
create_fasta(B_df, path + 'B_kmers.fasta')

A_df.to_csv(path + 'A_kmers.tsv', sep='\t')
B_df.to_csv(path + 'B_kmers.tsv', sep='\t')
print('created A_kmers.fasta and B_kmers.fasta')