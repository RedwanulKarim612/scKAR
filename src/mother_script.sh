#!/bin/bash

. ./config.env

ls $INPUT_DIR

dataset_path="$INPUT_DIR"

size=0
parallel_instances=0
declare -a pids=()
max_parallel=5

function check_pids () {
	for pid in "${pids[@]}"; do
		wait "$pid"
	done
	parallel_instances=0
	pids=()
}

function check_max_parallel () {
	if [ $parallel_instances -eq $max_parallel ]; then
		check_pids
	fi
}


start_global=`date +%s`


# files=$(find $dataset_path -maxdepth 1 -name "*.fastq.gz" | sort -g)

start=`date +%s`
bash jellyfish_count.sh $dataset_path/reads
end=`date +%s`
runtime=$((end-start))

echo "jellyfish count took: $runtime seconds"
echo "starting filtering"

start=`date +%s`
fas=$(find $dataset_path/jellyfish -name "*.fa" | sort -g)
for fa in $fas; do
    if [ ! -f "$fa.1_filtered.csv" ]; then
        echo "filtering $fa"
        bash faToCSV.sh $fa && python3 filter_kmers.py $fa.csv 31 && rm $fa $fa.csv &
        pids+=($!)
        parallel_instances=$((parallel_instances+1))
        check_max_parallel
    fi
done
end=`date +%s`
runtime=$((end-start))

echo "filtering 1-mers took: $runtime seconds"
echo "----------------------"

end_global=`date +%s`
runtime=$((end_global-start_global))
runtime_minutes=$((runtime/60))
echo "total runtime: $runtime_minutes minutes"
echo "total size: $size bytes"