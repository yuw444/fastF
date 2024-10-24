#!/bin/bash

module load umi_tools/1.1.2

N=0.5
R=0.5

# mkdir -p /scratch/g/chlin/Yu/FastQDesign/AIBM/benchmark/N_${N}_R_${R}
cd /scratch/g/chlin/Yu/FastQDesign/AIBM/benchmark/N_${N}_R_${R}

## step3 generate the feature matrix
umi_tools count --per-gene \
     --gene-tag=GX \
     --extract-umi-method=tag \
     --per-cell \
     --cell-tag=CB \
     --umi-tag=UB \
     -I filtered_genome_bam_N_${N}_R_${R}.bam \
     -S counts.tsv.gz
