#!/bin/sh

path=$1
files=$(find $path -name "*.fq.gz" -o -name "*.fastq.gz" | sort -g)

function kmer_count () {
    echo "------"
    echo "kmer_count of: $1"
    dir=$(dirname "$file")
    dir=$dir/../jellyfish
    base=$(basename "$file" .fq.gz)

    jf_path="${dir}/${base}.jf"
    fa_path="${dir}/${base}.fa"

    if [ -f $fa_path ]; then
        return
    fi

    if [ -f $fa_path.csv_1_filtered.csv ]; then
        return
    fi

    gunzip -c $1 | jellyfish count -m 31 -s 100M -t 8 -o $jf_path /dev/fd/0
    jellyfish dump $jf_path > $fa_path
    rm $jf_path
    echo "------"
}

i=0
for file in $files
do
    if [ ! -d $path/../jellyfish ]; then
        mkdir -v $path/../jellyfish
    fi
    if [[ $file == *"unmatched"* ]]; then
        continue
    fi
    kmer_count $file

    i=$((i+1))
done