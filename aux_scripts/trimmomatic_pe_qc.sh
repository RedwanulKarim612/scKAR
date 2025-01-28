#!/bin/sh

path=$1

files=$(find $path/fastqs/ -name "*_1.fastq.gz")


for file in $files; do
    filename=${file%_1.fastq.gz}
    if [ ! -f ${filename}_forward_paired.fq.gz ]; then
        echo "Trimming $filename"
        java -jar ./Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -phred33 \
        $file ${filename}_2.fastq.gz \
        ${filename}_forward_paired.fq.gz ${filename}_forward_unpaired.fq.gz \
        ${filename}_reverse_paired.fq.gz ${filename}_reverse_unpaired.fq.gz \
        ILLUMINACLIP:./Trimmomatic/adapters/NexteraPE-PE.fa:2:30:10:8:true \
        SLIDINGWINDOW:4:15 \
        MINLEN:31
    fi
done
