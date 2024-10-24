## Download the benchmark dataset from GEO

Accession "GSE269611"

## Using Cellranger to get BAM file and UMI matrix

```
cellranger count --id AI \
                 --fastqs /path/to/rawFastq \
                 --sample=AI \
                 --transcriptome=/path/to/library/mm10
```

## Benchmark

* **fastF**: fastF_disk.sh
* **Current**: 
  * step1_downsample.sh
  * step2_conversion.sh
