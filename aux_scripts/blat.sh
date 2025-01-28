#!/bin/sh


blat_path=$1
db_path=$2
contig_file=$3
output_file=$4
blast9_flag=$5

echo "running blat on $contig_file"

if [ "$blast9_flag" = "blast9" ] ; then
    echo "BLAST9 MODE"
    $blat_path $db_path $contig_file -out=blast9 $output_file
else 
    echo "BED MODE"
    $blat_path $db_path $contig_file $output_file
fi
