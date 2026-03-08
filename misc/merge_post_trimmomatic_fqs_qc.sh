#/bin/sh

path=$1

files=$(find $path/fastqs/ -name "*_1.fastq.gz")

for file in $files; do
    filename=${file%_1.fastq.gz}
    # get the file name removing its path
    basename=$(basename $filename)
    echo "Merging $filename"
    # Merge forward and reverse paired reads and save in another folder
    mkdir -p $path/merged_fastqs
    cat ${filename}_forward_paired.fq.gz ${filename}_reverse_paired.fq.gz > $path/merged_fastqs/${basename}.fq.gz
done