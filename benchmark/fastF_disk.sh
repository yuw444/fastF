#!/bin/bash

N=0.5
R=0.5

mkdir -p ./benchmark/N_${N}_R_${R}/fastF
cd ./benchmark/N_${N}_R_${R}/fastF

fastF bam2db -b /path/to/possorted_genome_bam.bam \
    -f /path/to/filtered_feature_bc_matrix/features.tsv.gz \
    -a /path/to/filtered_feature_bc_matrix/barcodes.tsv.gz \
    -d ./ai.db \
    -c 0.5 \
    -r 0.5 \
    -o .   \
    -s 926

