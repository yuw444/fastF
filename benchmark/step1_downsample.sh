#!/bin/bash

N=0.5
R=0.5

mkdir -p ./benchmark/N_${N}_R_${R}
cd ./benchmark/N_${N}_R_${R}

## step1 sample 50% of the cells
zcat /path/to/filtered_feature_bc_matrix/barcodes.tsv.gz | \
    awk -v n="$N" 'BEGIN {srand(926)} {if (rand() <= n) print}' > cells_$N.txt


## step2 sample 50% of the reads
module load samtools/1.20

awk 'NR==FNR {barcodes[$1]; next} 
    $0 ~ /^@/ || 
    ($0 ~ /CB:Z:/ && $0 ~ /UB:Z:/ && ($0 ~ /xf:i:25/ || $0 ~ /xf:i:17/) && 
        substr($0, index($0, "CB:Z:") + 5, 18) in barcodes)' \
    cells_${N}.txt \
    <(samtools view -@ 1 -h -s $R /path/to/possorted_genome_bam.bam) | \
    samtools view -@ 1 -bo ./filtered_genome_bam_N_${N}_R_${R}.bam

samtools index filtered_genome_bam_N_${N}_R_${R}.bam
