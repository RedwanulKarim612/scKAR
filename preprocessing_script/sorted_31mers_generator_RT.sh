#!/bin/sh

# Argument takes the genome.fa.gz file and the output file path

gunzip -c $1 | jellyfish count -m 31 -s 100M -t 8 -o $1.jf /dev/fd/0
jellyfish dump $1.jf > $1.fa

awk '!/^>/' $1.fa | sort > $2
