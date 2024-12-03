import pandas as pd
import sys
import os

def cal_tpm(df, kmer_len, total_kmer_count):
    df['TPM'] = (df['count'] / total_kmer_count) * 10**6
    return df['TPM']

def filter_1_count(file):
    df = pd.read_csv(file, names=['kmer', 'count'], header=None)
    df = df[df['count'] > 1]
    return df


filename = sys.argv[1]
dirname  = os.path.dirname(filename)

filtered_df = filter_1_count(filename)
total_kmer_count = filtered_df['count'].sum()

tpm = cal_tpm(filtered_df, int(sys.argv[2]), total_kmer_count)
norm_df = pd.DataFrame({'kmer': filtered_df['kmer'], 'tpm': tpm})

filename = filename.replace('.fastq.gz.fa.csv', '')
norm_df.to_csv(filename+"_1_filtered.csv", index=False, header=False)

with open(dirname + '/../tpm_sum.csv', 'a') as f:
    only_filename = os.path.basename(filename)
    f.write(only_filename + '\t' + str(total_kmer_count)+'\n')
    f.flush()

