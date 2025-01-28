import pandas as pd
import gffutils
import sys
# path = '../../data/axolotl/1_wpi_demux/final_results/bipartition_0/abyss/'
# db_path = '/media/supertramp/ssd_ratul/thesis/reference_genome_axolotl/axolotl_omics/AmexT_v47_v6DD.db'
path = sys.argv[1]
db_path = sys.argv[2]
sam_df = pd.read_csv(path + 'B_contigs_alignment.sam', header=None, sep='\t', comment='@', usecols=[0,1,2,3,4,5,6,7,8,9])   
sam_df.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ']
sam_df['contig_length'] = sam_df['SEQ'].str.len()
sam_df.sort_values(by=['contig_length'], inplace=True, ascending=False)      
# sam_df.head()
# bed_df = pd.read_csv('/media/supertramp/supertramp_s/Thesis/data/axolotl/1_wpi_demux/final_results/sand/A_contigs_alignment.bed', header=None, sep='\t')   
# bed_df.columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
# bed_df.sort_values(by=['name'], inplace=True, ascending=False)
sam_df.sort_values(by=['QNAME'], inplace=True, ascending=False)
# bed_df.head()
# sam_df.head()
def calculate_end_position(start_pos, cigar_string):
    # Operations that consume the reference sequence
    consume_reference_ops = set('MDN=X')
    length_consumed = 0
    num_buffer = ''  # Buffer to accumulate digits

    for char in cigar_string:
        if char.isdigit():
            num_buffer += char  # Build up the number
        elif char in consume_reference_ops:
            length_consumed += int(num_buffer)  # Add the length for operations consuming the reference
            num_buffer = ''  # Reset the number buffer
        else:
            # For operations that do not consume the reference sequence, just reset the buffer
            num_buffer = ''

    # Calculate the end position
    end_position = start_pos + length_consumed - 1
    return end_position

# Example usage
# start_position_example = 100
# cigar_string_example = "10M2I5M3D7M"
# calculate_end_position(start_position_example, cigar_string_example)
sam_df['END'] = sam_df.apply(lambda row: calculate_end_position(row['POS'], row['CIGAR']), axis=1)
# sam_df.head()
# db_path = '/media/supertramp/ssd_ratul/thesis/reference_genome_axolotl/axolotl_omics/'

# db = gffutils.FeatureDB(db_path + "AmexT_v47_v6DD.db")
db = gffutils.FeatureDB(db_path)
replace_map= {'chr1p':'CM010927.2',
'chr1q':'CM010928.2',
'chr2p':'CM010929.2',
'chr2q':'CM010930.2',
'chr3p':'CM010931.2',
'chr3q':'CM010932.2',
'chr4p':'CM010933.2',
'chr4q':'CM010934.2',
'chr5p':'CM010935.2',
'chr5q':'CM010936.2',
'chr6p':'CM010937.2',
'chr6q':'CM010938.2',
'chr7p':'CM030316.1',
'chr7q':'CM030317.1',
'chr8p':'CM030318.1',
'chr8q':'CM030319.1',
'chr9p':'CM030320.1',
'chr9q':'CM030321.1',
'chr10p':'CM030322.1',
'chr10q':'CM030323.1',
'chr11p':'CM030324.1',
'chr11q':'CM030325.1',
'chr12p':'CM030326.1',
'chr12q':'CM030327.1',
'chr13p':'CM030328.1',
'chr13q':'CM030329.1',
'chr14p':'CM030330.1',
'chr14q':'CM030331.1'
}
# reverse the dictionary
replace_map = {v: k for k, v in replace_map.items()}
sam_df['RNAME'] = sam_df['RNAME'].replace(replace_map)
sam_df
annotation_map = {
    'contig_id' :[],
    'contig': [],
    'chromosome' :[],
    'start' :[],
    'end' :[],
    'gene_id' :[],
    'gene_name' :[],
    'transcript_id' :[],
    'flag':[],
    'mapq' :[],
    'cigar' :[],
    'feature_type' :[],
    'homolog': [],
    'orf_type':[],              
    'CDS': [],
    'exon_number':[],
    'peptide':[],                                        
}
def insert_into_map(contig_id, contig, flag, chromosome='*', start=-1, end=-1, gene_id='*', gene_name='*', transcript_id='*', mapq=0, cigar='*', feature_type='*', homolog='*', orf_type='*', CDS='*', exon_number=-1, peptide='*'):
    annotation_map['contig_id'].append(contig_id)
    annotation_map['contig'].append(contig)
    annotation_map['chromosome'].append(chromosome)
    annotation_map['start'].append(start)
    annotation_map['end'].append(end)
    annotation_map['gene_id'].append(gene_id)
    annotation_map['gene_name'].append(gene_name)
    annotation_map['transcript_id'].append(transcript_id)                                                                                       
    annotation_map['flag'].append(flag)
    annotation_map['mapq'].append(mapq)
    annotation_map['cigar'].append(cigar)
    annotation_map['feature_type'].append(feature_type)
    annotation_map['homolog'].append(homolog)
    annotation_map['orf_type'].append(orf_type)    
    annotation_map['CDS'].append(CDS)
    annotation_map['exon_number'].append(exon_number)
    annotation_map['peptide'].append(peptide)

# f_types = set()
for index, row in sam_df.iterrows():
    contig_id = row['QNAME']
    start = row['POS']
    end = row['END']
    chrom = row['RNAME']
    #check if contig is unmapped
    if chrom == '*':
        insert_into_map(contig_id, row['SEQ'], flag=row['FLAG'])
        continue
    # print(contig, start, end, chrom)
    for feature in db.region(region=(chrom, start, end)):
                                                    
        transcript_id = '*'
        gene_id = '*'                                               
        gene_name = '*'
        homolog = '*'
        orf_type = '*'
        cds = '*'
        exon_number = -1
        peptide = '*'


        if 'gene_id' in feature.attributes:
            gene_id = feature.attributes['gene_id'][0]
        if 'gene_name' in feature.attributes:
            gene_name = feature.attributes['gene_name'][0]
        if 'transcript_id' in feature.attributes:
            transcript_id = feature.attributes['transcript_id'][0]
        if 'homolog' in feature.attributes:
            homolog = feature.attributes['homolog'][0]
        if 'ORF_type' in feature.attributes:
            orf_type = feature.attributes['ORF_type'][0]                                                
        if 'CDS' in feature.attributes:
            cds = feature.attributes['CDS'][0]
        if 'exon_number' in feature.attributes:
            exon_number = feature.attributes['exon_number'][0]
        if 'peptide' in feature.attributes:
            peptide = feature.attributes['peptide'][0]                                                      
        insert_into_map(contig_id, row['SEQ'], flag=row['FLAG'], chromosome=chrom, start=start, end=end, 
                            gene_id=gene_id, 
                            gene_name=gene_name, 
                            transcript_id=transcript_id,
                            mapq=row['MAPQ'], 
                            cigar=row['CIGAR'], 
                            feature_type=feature.featuretype, 
                            homolog=homolog,
                            orf_type=orf_type,
                            CDS=cds,
                            exon_number=exon_number,
                            peptide=peptide
                         )
        # break

    # break
# print(f_types)
annotation_df = pd.DataFrame(annotation_map)
annotation_df.head()
annotation_df['homolog'].unique()


annotation_df.to_csv(path + 'B_contigs_annotation.tsv', index=False, sep='\t')







