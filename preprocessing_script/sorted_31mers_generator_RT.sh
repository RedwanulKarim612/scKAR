#!/bin/sh

# Argument takes the genome.fa file.

gunzip -c $1 | jellyfish count -m 31 -s 100M -t 8 -o $jf_path /dev/fd/0
jellyfish dump $jf_path > $fa_path

awk '!/^>/' $1 | sort > 31_mers_only.txt
