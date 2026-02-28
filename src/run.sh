#!/bin/sh

. ./config.env

# k-mer abundance generation
cd ./kmer-filtering
bash ./mother_script.sh $INPUT_DIR $FILTER_REFERENCE_KMERS $REFERENCE_TRANSCRIPTOME_SORTED_31MERS_PATH $NUMBER_OF_THREADS

# condition generation
if [[ "$MODE" == "GENE_EXPRESSION" ]] then
    cd ../clustering
    python3 graph-clustering.py "$INPUT_DIR"/expression_matrix/expression.csv $CLUSTERING_ALGO $MIN_GENES $MIN_CELLS $N_NEIGHBORS $N_PCS $RESOLUTION
else if [[ "$MODE" == "CUSTOM_METADATA" ]] then
    cd ../clustering
    python3 generate_bipartitions_from_metadata.py "$INPUT_DIR"
else if [[ "$MODE" == "KMER_ABUNDANCE" ]] then
    total_cell_count=$(wc -l "$INPUT_DIR"/tpm_sum.csv | awk '{print $1}')
    cd ../clustering
    g++ merge_adj_for_clustering.cpp -o merge_adj_for_clustering
    ./merge_adj_for_clustering $INPUT_DIR/adj $INPUT_DIR/expression_matrix/kmer_abundance.csv $total_cell_count
    python3 graph-clustering.py "$INPUT_DIR"/expression_matrix/kmer_abundance.csv $CLUSTERING_ALGO $MIN_GENES $MIN_CELLS $N_NEIGHBORS $N_PCS $RESOLUTION
else
    echo "Invalid mode"
    exit 1
fi
# F-test
cd ../f-test
python3 preprocess_clustering_results.py "$INPUT_DIR"/tpm_sum.csv "$INPUT_DIR"/clustering_results/cluster.csv  "$INPUT_DIR"/cluster_tpm.csv
g++ adj_to_mat.cpp -o f_test.out
./f_test.out $INPUT_DIR $PSEUDOBULK_SIZE $MIN_ROW_THREDSHOLD

# DE test
cd ../de-test
bash ./deseq_runner.sh $INPUT_DIR $MIN_ROW_COUNT $MIN_COL_COUNT $NUMBER_OF_THREADS $LOG2FC $PVAL $BASE_MEAN_THRESHOLD

# K-mer Assembly
cd ../merge_kmers
bipartition_folders=$(find $INPUT_DIR/final_results -mindepth 1 -maxdepth 1 -type d)
for folder in $bipartition_folders; do
    bash ./create_contigs.sh $folder $BLAT_PATH "$INPUT_DIR"/genome_assembly/genome.2bit
done
