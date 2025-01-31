#!/bin/sh
path=$1
blat_path=$2
db_path=$3
log2fc=$4
pval=$5
base_mean_threshold=$6
# if blast9_flag is true, then use blast9 flag in blat
# if blast9_flag is false, then don't use blast9 flag in blat
echo "blast9 flag is true"
python3 ./merge_deseq_results.py $path $log2fc $pval $base_mean_threshold
# # create a folder named dekupl to store the output of dekupl if it doesn't exist
# if [ ! -d $path/dekupl ]; then
#     mkdir -p $path/dekupl
# fi
# if [ ! -d $path/spades ]; then
#     mkdir -p $path/spades
# fi
if [ ! -d $path/abyss ]; then
    mkdir -p $path/abyss
fi
# echo "-------------Merging kmers using dekupl----------------"
# ./dekupl_merge_tags/mergeTags -k 31 -m 19 -n $path/filtered_kmers.tsv > $path/dekupl/dekupl_contigs.tsv
# python3 ./create_fastas_dekupl.py $path/dekupl/
# ./blat.sh $blat_path $db_path $path/dekupl/A_contigs.fasta $path/dekupl/A_contigs_alignment.psl "blast9"
# ./blat.sh $blat_path $db_path $path/dekupl/B_contigs.fasta $path/dekupl/B_contigs_alignment.psl "blast9"
# echo "-------------Finished merging kmers using dekupl----------------"

# echo "-------------Merging kmers using spades----------------"
# >> spades.log
# spades --sc --only-assembler -t 6 -k 19,21,23 -s $path/A_kmers.fasta -o $path/spades/assembled_sc_k19_k21_k23_A >> spades.log
# spades --sc --only-assembler -t 6 -k 19,21,23 -s $path/B_kmers.fasta -o $path/spades/assembled_sc_k19_k21_k23_B >> spades.log
# python3 ./create_fastas_spades.py $path/spades/
# ./blat.sh $blat_path $db_path $path/spades/A_contigs.fasta $path/spades/A_contigs_alignment.psl "blast9"
# ./blat.sh $blat_path $db_path $path/spades/B_contigs.fasta $path/spades/B_contigs_alignment.psl "blast9"
# echo "-------------Finished merging kmers using spades----------------"

echo "-------------Merging kmers using abyss----------------"
>> abyss.log
abyss-pe k=25 name=A se=$path/A_kmers.fasta >> abyss.log
abyss-pe k=25 name=B se=$path/B_kmers.fasta >> abyss.log
cat A-unitigs.fa > $path/abyss/A_contigs.fasta
rm A-*
cat B-unitigs.fa > $path/abyss/B_contigs.fasta
rm B-*
python3 create_tsv_from_fasta.py $path/abyss/

./blat.sh $blat_path $db_path $path/abyss/A_contigs.fasta $path/abyss/A_contigs_alignment.psl "blast9"
./blat.sh $blat_path $db_path $path/abyss/B_contigs.fasta $path/abyss/B_contigs_alignment.psl "blast9"