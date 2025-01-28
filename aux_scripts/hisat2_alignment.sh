#!/bin/sh

fa_path=$1
index_path='../../../../ssd_ratul/thesis/axolotl_ref_T_index/hisat2_index_refT_axolotl'
output_path=$2

hisat2 -x $index_path -f $fa_path -S $output_path
