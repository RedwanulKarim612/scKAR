#!/bin/bash

input_file=$1

awk 'NR%2==1{sub(/^>/, "", $0); count = $0; next} {print $0 "," count}' $input_file > $input_file.csv