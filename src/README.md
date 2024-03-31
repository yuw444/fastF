

# ds bam 

* database in disk

```
./fastF bam2db -b /scratch/g/chlin/Yu/fastF/data/test.bam \
    -f /scratch/g/chlin/Yu/fastF/data/features.tsv.gz \
    -a /scratch/g/chlin/Yu/fastF/data/barcodes.txt.gz \
    -d sample.db \
    -c 0.5 \
    -r 0.5 \
    -o . \
    -s 926 

```
* database in memory
```
./fastF bam2db -b /scratch/g/chlin/Yu/fastF/data/test.bam \
    -f /scratch/g/chlin/Yu/fastF/data/features.tsv.gz \
    -a /scratch/g/chlin/Yu/fastF/data/barcodes.txt.gz \
    -d :memory: \
    -c 1.0 \
    -r 1.0 \
    -o . \
    -s 926 

```