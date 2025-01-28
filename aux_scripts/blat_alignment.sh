#!/bin/sh

blat_binary_path=$1
contigs_path=$2
output_path=$3

$1/blat $1/hg38.2bit $2 -out=blast9 $3/contigs_vs_hg38.blast9.psl