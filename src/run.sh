#!/bin/sh

. ./config.env

cd ./kmer-filtering
bash ./mother_script.sh $INPUT_DIR

cd ..
cd ./f-test
g++ adj_to_mat.cpp -o f_test.out
./f_test.out $INPUT_DIR 70000 100

cd ..
cd ./de-test
bash ./deseq_runner.sh $INPUT_DIR