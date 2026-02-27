#!/bin/sh

dataset_path=$1
ref_filter=$2
ref_transcriptome_file=$3
num_threads=$4

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

start=`date +%s`
bash ./jellyfish_count.sh $dataset_path/reads
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

echo "creating adjacency list"
start=`date +%s`

g++ create_adj_list.cpp -o create_adj_list
./create_adj_list $dataset_path

end=`date +%s`
runtime=$((end-start))
echo "creating adjacency list took: $runtime seconds"

if [[$ref_filter == true]]; then
	echo "filtering references"
	start=`date +%s`

	g++ ref_filter.cpp -o ref_filter
	./ref_filter $dataset_path $ref_transcriptome_file $num_threads
	end=`date +%s`

	runtime=$((end-start))
	echo "filtering references took: $runtime seconds"