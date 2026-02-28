#!/bin/sh

path=$1
merge_method=$2

if [ "$merge_method" = "hmp" ]; then
  echo "Using Harmonic Mean P-value (HMP) method for merging p-values."
  Rscript hmp.R "$path/abyss/A_contigs.tsv" "$path/abyss/A_kmers.tsv" "$path/abyss/A_contigs_pvalues.tsv" 
  Rscript hmp.R "$path/abyss/B_contigs.tsv" "$path/abyss/B_kmers.tsv" "$path/abyss/B_contigs_pvalues.tsv" 
  
elif [ "$merge_method" = "fisher" ]; then
  python fisher.py "$path/abyss/A_contigs.tsv" "$path/abyss/A_kmers.tsv" "$path/abyss/A_contigs_pvalues.tsv" "$merge_method"
  python fisher.py "$path/abyss/B_contigs.tsv" "$path/abyss/B_kmers.tsv" "$path/abyss/B_contigs_pvalues.tsv" "$merge_method"
elif [ "$merge_method" = "stouffer" ]; then
  python fisher.py "$path/abyss/A_contigs.tsv" "$path/abyss/A_kmers.tsv" "$path/abyss/A_contigs_pvalues.tsv" "$merge_method"
  python fisher.py "$path/abyss/B_contigs.tsv" "$path/abyss/B_kmers.tsv" "$path/abyss/B_contigs_pvalues.tsv" "$merge_method"
else
  echo "Unsupported merge method: $merge_method. Please choose 'hmp', 'fisher', or 'stouffer'."
  exit 1
fi