#!/bin/sh

path=$1
db_path=$2
cut_len=$3

echo "Annotating contigs in $path/A_contigs_alignment.psl"
python3 annotation_gffutils.py $db_path $path/A_contigs_alignment.psl $path/A_contigs.tsv $cut_len
# echo "Annotating contigs in $path/B_contigs_alignment.psl"
# python3 annotation_gffutils.py $db_path $path/B_contigs_alignment.psl $path/B_contigs.tsv $cut_len
