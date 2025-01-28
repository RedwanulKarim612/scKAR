#!/bin/sh

file=$1
folder=$(dirname $file)
filename=$(basename $file)
contig=${filename//kmers.fasta/contigs_abyss.fasta}

if [ ! -d $folder/${contig}_abyss_k25 ]; then
    mkdir $folder/${contig}_abyss_k25
fi

cp $file $folder/${contig}_abyss_k25/

k

# mv $folder/${contig}_abyss_k25/abyss-unitigs.fa $folder/$contig
