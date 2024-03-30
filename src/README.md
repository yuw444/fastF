

# ds bam 

* database in disk

```
./fastF bam2db -b /scratch/u/yu89975/fastF/data/test.bam \
    -f /scratch/u/yu89975/fastF/data/features.tsv.gz \
    -a /scratch/u/yu89975/fastF/data/barcodes.txt.gz \
    -d sample.db \
    -c 1.0 \
    -r 1.0 \
    -o . \
    -s 926 

```
* database in memory
```
./fastF bam2db -b /scratch/u/yu89975/fastF/data/test.bam \
    -f /scratch/u/yu89975/fastF/data/features.tsv.gz \
    -a /scratch/u/yu89975/fastF/data/barcodes.txt.gz \
    -d :memory: \
    -c 1.0 \
    -r 1.0 \
    -o . \
    -s 926 

```