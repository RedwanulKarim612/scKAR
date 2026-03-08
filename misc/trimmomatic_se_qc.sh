#!/bin/sh

path=$1

# Find all fastq files that have not been previously trimmed
files=$(find $path/fastqs/ -name "*.fastq.gz")

TRIMMOMATIC_JAR="../../Tools/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar"

ADAPTER_FILE="../../Tools/Trimmomatic/adapters/chromium_10x_3prime_v3.fa"

# makedir named qc_fastqs
if [ ! -d "$path/qc_fastqs" ]; then
    mkdir "$path/qc_fastqs"
fi 

# Loop over the found files
for file in $files; do
    filename=$(basename "$file" .fq.gz)
    trimmed_file="${path}/qc_fastqs/${filename}_trimmed.fq.gz"

    if [ ! -f "$trimmed_file" ]; then
        echo "Trimming $file"
        java -jar "$TRIMMOMATIC_JAR" SE -threads 4 -phred33 \
        "$file" \
        "$trimmed_file" \
        ILLUMINACLIP:${ADAPTER_FILE}:2:30:20 \
        SLIDINGWINDOW:4:20 \
        MINLEN:31 # Minimum length of reads to keep after trimming
        # Seed mismatches: 2, Simple clip threshold: 15
        # Window size: 4, Required quality: 10
    else
        echo "Trimmed file already exists: $trimmed_file"
    fi
done
