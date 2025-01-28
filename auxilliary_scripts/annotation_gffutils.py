import pandas as pd
import gffutils
import sys

db_path = sys.argv[1]
psl_file = sys.argv[2]
contigs_file = sys.argv[3]
cut_length = int(sys.argv[4])

db = gffutils.FeatureDB(db_path + "GENCODE_v45.db")
polyA_db = gffutils.FeatureDB(db_path + "polyA_v45.db")
tRNA_db = gffutils.FeatureDB(db_path + "tRNA_v45.db")
psl_df = pd.read_csv(psl_file, delimiter='\t', comment='#', header=None, names=['name', 'subject', 'identity', 'length', 'mismatches', 'gap_openings', 'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score'])
contigs_df = pd.read_csv(contigs_file, sep='\t')
print('contigs count: ', contigs_df.head)
# plt.hist(contigs_df['length'], bins=[i for i in range(40, 50, 1)])
contigs_df['contig_length'] = contigs_df['contig'].apply(lambda x: len(x))
contigs_df = contigs_df[contigs_df['contig_length'] >= cut_length]
psl_df = pd.merge(psl_df, contigs_df, on='name')
# psl_df.head
psl_df = psl_df[psl_df['contig_length']==psl_df['length']]
# psl_df.shape
psl_df.sort_values(by='e_value', inplace=True)
psl_df.drop_duplicates(subset=['subject', 's_start', 's_end'], inplace=True, keep='first')
print('psl count: ', psl_df.shape[0]    )
lncRNAs = {
    'name': [],
    'contig' : [],
    'chr': [],
    'start': [],
    'end': [],
    'gene_name': [],
}
pseudoGenes = {
    'name': [],
    'contig' : [],
    'chr': [],
    'start': [],
    'end': [],
    'gene_name': [],
}
retained_introns = {
    'name': [],
    'contig' : [],
    'chr': [],
    'start': [],
    'end': [],
    'gene_name': [],
}
transcribed_unprocessed_pseudogenes = {
    'name': [],
    'contig' : [],
    'chr': [],
    'start': [],
    'end': [],
    'gene_name': [],
}
miRNAs = {
    'name': [],
    'contig' : [],
    'chr': [],
    'start': [],
    'end': [],
    'gene_name': [],
}
polyAs = {
    'name': [],
    'chr': [],
    'start': [],
    'end': [],
    'gene_name': [],
}
tRNAs = {
    'name': [],
    'contig' : [],
    'chr': [],
    'start': [],
    'end': [],
    'gene_name': [],
}
def add_to_annotation_dict(annotation_dict, feature, name, contig chr, start, end):
    annotation_dict['name'].append(name)
    annotation_dict['contig'].append(contig)
    annotation_dict['chr'].append(chr)
    annotation_dict['start'].append(start)
    annotation_dict['end'].append(end)
    # annotation_dict['strand'].append(feature.strand)
    annotation_dict['gene_name'].append(feature.attributes['gene_name'][0])

cnt = 0
for index, row in psl_df.iterrows():
    chr, start, end = row['subject'], row['s_start'], row['s_end']
    features = db.region(region=(chr, start, end))
    for f in features:
        # print(f)
        if 'transcript_type' in f.attributes:
            if f.attributes['transcript_type'][0]=='retained_intron':
                add_to_annotation_dict(retained_introns, f, row['name'], row['contig'], chr, start, end)
        if f.attributes['gene_type'][0]=='lncRNA':
            add_to_annotation_dict(lncRNAs, f, row['name'], row['contig'], chr, start, end)
        if f.attributes['gene_type'][0]=='processed_pseudogene':
            add_to_annotation_dict(pseudoGenes, f, row['name'], row['contig'], chr, start, end)
        if f.attributes['gene_type'][0]=='transcribed_processed_pseudogene':
            add_to_annotation_dict(transcribed_unprocessed_pseudogenes, f, row['name'], row['contig'], chr, start, end)
        if f.attributes['gene_type'][0]=='miRNA':
            add_to_annotation_dict(miRNAs, f, row['name'], row['contig'], chr, start, end)
    cnt+=1
    if cnt%1000==0:
        print(cnt)
for index, row in psl_df.iterrows():
    chr, start, end = row['subject'], row['s_start'], row['s_end']
    features = polyA_db.region(region=(chr, start, end))
    for f in features:
        add_to_annotation_dict(polyAs, f, row['name'], row['contig'], chr, start, end)
    
for index, row in psl_df.iterrows():
    chr, start, end = row['subject'], row['s_start'], row['s_end']
    features = tRNA_db.region(region=(chr, start, end))
    for f in features:
        add_to_annotation_dict(tRNAs, f, row['name'], row['contig'], chr, start, end)
    
lncRNAs_df = pd.DataFrame(lncRNAs)
# lncRNAs_df.drop_duplicates(subset=['name'], inplace=True, keep='first')
# lncRNAs_df.drop_duplicates(subset=['gene_name'], inplace=True, keep='first')
lncRNAs_df.to_csv(contigs_file.replace('.tsv', '_lncRNA.tsv'), sep='\t', index=False) 

retained_introns_df = pd.DataFrame(retained_introns)
# retained_introns_df.drop_duplicates(subset=['name'], inplace=True, keep='first')
# retained_introns_df.drop_duplicates(subset=['gene_name'], inplace=True, keep='first')
retained_introns_df.to_csv(contigs_file.replace('.tsv', '_retained_introns.tsv'), sep='\t', index=False)

pseudoGenes_df = pd.DataFrame(pseudoGenes)
# pseudoGenes_df.drop_duplicates(subset=['name'], inplace=True, keep='first')
# pseudoGenes_df.drop_duplicates(subset=['gene_name'], inplace=True, keep='first')
pseudoGenes_df.to_csv(contigs_file.replace('.tsv', '_pseudoGenes.tsv'), sep='\t', index=False)

transcribed_unprocessed_pseudogenes_df = pd.DataFrame(transcribed_unprocessed_pseudogenes)
transcribed_unprocessed_pseudogenes_df.drop_duplicates(subset=['name'], inplace=True, keep='first')
transcribed_unprocessed_pseudogenes_df.drop_duplicates(subset=['gene_name'], inplace=True, keep='first')
transcribed_unprocessed_pseudogenes_df.to_csv(contigs_file.replace('.tsv', '_transcribed_unprocessed_pseudogenes.tsv'), sep='\t', index=False)

polyAs_df = pd.DataFrame(polyAs)
polyAs_df.drop_duplicates(subset=['name'], inplace=True, keep='first')
polyAs_df.drop_duplicates(subset=['gene_name'], inplace=True, keep='first')
polyAs_df.to_csv(contigs_file.replace('.tsv', '_polyAs.tsv'), sep='\t', index=False)

tRNAs_df = pd.DataFrame(tRNAs)
tRNAs_df.drop_duplicates(subset=['name'], inplace=True, keep='first')
tRNAs_df.drop_duplicates(subset=['gene_name'], inplace=True, keep='first')
tRNAs_df.to_csv(contigs_file.replace('.tsv', '_tRNAs.tsv'), sep='\t', index=False)
# find intersections between lncRNAs and retained introns

miRNAs_df = pd.DataFrame(miRNAs)
miRNAs_df.drop_duplicates(subset=['name'], inplace=True, keep='first')
miRNAs_df.drop_duplicates(subset=['gene_name'], inplace=True, keep='first')
miRNAs_df.to_csv(contigs_file.replace('.tsv', '_miRNAs.tsv'), sep='\t', index=False)