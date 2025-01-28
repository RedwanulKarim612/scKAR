#!/bin/sh

. ./config.env

# k-mer abundance generation
cd ./kmer-filtering
bash ./mother_script.sh $INPUT_DIR

# condition generation
cd ../clustering
python3 graph-clustering.py "$INPUT_DIR"/expression_matrix/expression.csv $CLUSTERING_ALGO $MIN_GENES $MIN_CELLS $N_NEIGHBORS $N_PCS $RESOLUTION

# F-test
cd ../f-test
python3 preprocess_clustering_results.py "$INPUT_DIR"/tpm_sum.csv "$INPUT_DIR"/clustering_results/cluster.csv  "$INPUT_DIR"/cluster_tpm.csv
g++ adj_to_mat.cpp -o f_test.out
./f_test.out $INPUT_DIR 70000 100

# DE test
cd ../de-test
bash ./deseq_runner.sh $INPUT_DIR

# K-mer Assembly
