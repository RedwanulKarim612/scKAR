#!/bin/sh

. ./config.env

cd ./kmer-filtering
bash ./mother_script.sh $INPUT_DIR

