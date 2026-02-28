import pandas as pd 
import scipy.stats as stats
import sys

contigs_file = sys.argv[1]
kmers_file   = sys.argv[2]
output_file  = sys.argv[3]
merge_method = sys.argv[4]

k = 31
EPS = 1e-200

def get_all_kmers(contig, k=31):
    kmers = []
    for i in range(len(contig) - k + 1):
        kmers.append(contig[i:i+k])
    return kmers    

def get_kmer_pvalues(kmers, kmers_df):
    pvalues = []
    for kmer in kmers:
        if kmer in kmers_df.index:
            p = kmers_df.loc[kmer].padj
            if pd.notna(p):
                pvalues.append(min(max(float(p), EPS), 1.0))

    if len(pvalues) == 0:
        return 1.0

    return stats.combine_pvalues(pvalues, method=merge_method)[1]

contigs_df = pd.read_csv(contigs_file, sep='\t')

contigs_df['kmers'] = contigs_df['contig'].apply(lambda x: get_all_kmers(x, k))

kmers_df = pd.read_csv(kmers_file, sep='\t', index_col=0)

kmers_df.index.names = ['kmer']

contigs_df['pvalue'] = contigs_df['kmers'].apply(lambda x: get_kmer_pvalues(x, kmers_df))

contigs_df = contigs_df[contigs_df['pvalue'] < 0.05] 
contigs_df = contigs_df.sort_values('pvalue')

contigs_df.to_csv(output_file, sep='\t', index=False, columns=['contig', 'pvalue'])