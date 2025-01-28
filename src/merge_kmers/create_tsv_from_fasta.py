import pandas as pd
from Bio import SeqIO
import sys

def fasta_to_tsv(input_fasta, output_tsv):
    data = []
    
    # Read FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequence_id = record.id
        sequence = str(record.seq)
        data.append([sequence_id, sequence])

    # Convert to DataFrame
    df = pd.DataFrame(data, columns=['name', 'contig'])
    
    # Write to TSV file
    df.to_csv(output_tsv, sep='\t', index=False)

fasta_to_tsv(sys.argv[1] + 'A_contigs.fasta', sys.argv[1] + 'A_contigs.tsv')
# fasta_to_tsv(sys.argv[1] + 'B_contigs.fasta', sys.argv[1] + 'B_contigs.tsv')