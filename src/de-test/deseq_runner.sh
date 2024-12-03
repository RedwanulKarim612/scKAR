#!/bin/sh

path=$1

metadata_folder="$path/bipartitions/"
matrix_folder="$path/f_test_results/"
tpm_file="$path/cluster_tpm.csv"
rowThreshold=$2

metadata=$(ls "$metadata_folder" | sort -g)
matrix=$(ls "$matrix_folder" | sort -g)

if [ ! -d "$path/deseq_results" ]; then
    mkdir "$path/deseq_results"
fi


for meta in $metadata; do
    for mat in $matrix; do
        # echo " --------------- "$meta"_"$mat" -----------------"
        if [ ! -f "$path/deseq_results/"$meta"_"$mat"_gp.csv" ]; then
            echo "Metadata file "$meta"_"$mat".csv does not exists"
            echo "Running DESeq2 for $meta and $mat"
            Rscript ./deseq2_glmgamPoi.R "$metadata_folder/$meta" "$matrix_folder/$mat" "$tpm_file" "$path"/deseq_results/"$meta"_"$mat"_gp.csv
            # break
        fi        
    done
    # break
done
