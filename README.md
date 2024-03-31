## fastF

`fastF` aims to **sample**, **filter**, and **summarise** single-cell RNA sequencing result, [FastQ](https://en.wikipedia.org/wiki/FASTQ_format) file. It consists 5 subcommands for the specific tasks.

 ```
> ./fastF --help
Usage: fastF <subcommands> [args]

subcommands:

        freq:    Find all the cell barcode whitelist and their frequencies.
        filter:  Filter fastq file using cell barcode whitelist and read depth.
        crb:     Extract CR and CB tags from bam file and summarize them with 
                     frequencies to a tsv file.
        bam2db:  Filter bam file with desired cell proportion and read depth, 
                     then summarise it UMI matrix.
        extract: Extract the tag of bam file and write it to tag_summary.csv.


    -h, --help
 ```   
    
### `freq`: Find all the cell barcode whitelist and their frequencies.

* Compare to [umi-tools whitelist](https://umi-tools.readthedocs.io/en/latest/reference/whitelist.html#usage), it extracts all the unique cell barcodes and their frequencies.
  
```
> ./fastF freq --help

Find all the cell barcode whitelist and their frequencies.

    -h, --help        

Basic options
    -R, --R1=<str>    path to R1 fastq files
    -o, --out=<str>   path to output whitelist
    -l, --len=<int>   length of cell barcode
    -u, --umi=<int>   length of UMI
```

### `filter`: Filter fastq file using cell barcode `whitelist` and `read depth`.

* It is the motivation of the entire project, see the similar [questions](https://www.biostars.org/post/search/?query=sample+cells+from+fastQ+file) arised in the community, and [FastQDesign](https://github.com/yuw444/FastQDesign). 

* When filter through the assumed `read depth`, it is compare to [seqtk](https://github.com/lh3/seqtk), but it could process all `I1`, `R1` and `R2` concurrently as needed. 

* Even more, it could filter FastQ reads according to the given cell barcode whitelist.

* Of course, `whitelist` and `read depth` can be set at the same time.

```
> ./fastF filter --help
Filter fastq file using cell barcode whitelist and read depth.

    -h, --help                

Basic options
    -I, --I1=<str>            optional, path to sample I1 fastq files
    -R, --R1=<str>            required, path to sample R1 fastq files
    -r, --R2=<str>            optional, path to sample R2 fastq files
    -o, --out=<str>           dir to output fastq files
    -w, --whitelist=<str>     whitelist of cell barcodes
    -l, --len=<int>           length of cell barcode
    -s, --seed=<int>          seed for random number generator
    -t, --rate=<flt>          rate of reads to keep after matching cell barcodes
    -a, --allcells            keep all reads with cell barcode
```

### `crb`:     Extract CR and CB tags from bam file and summarize them with frequencies to a tsv file.

* Once again, this command is very compared to [umi-tools whitelist](https://umi-tools.readthedocs.io/en/latest/reference/whitelist.html#usage). Instead of calculating the distance between all the raw cell barcodes, and clustering, it directly uses the tags `CR` and `CB` information in the [bam](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam) file produced by [cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest). 

* The result is summarise by the corrected cell barcode with its variants, along with their frequencies. 

```
> ./fastF crb --help

Extract CR and CB tags from bam file and summarize them with frequencies to a tsv file.

    -h, --help        
    -b, --bam=<str>   path to bam file
    -o, --out=<str>   path to output directory
```

### `bam2db`:  Filter bam file with desired cell proportion and read depth, then summarise it UMI matrix.

* `bam2db` is an upgraded version of `filter`, it can automatically find `whitelist` and subsample cell and read depth. 

* Compared to `filter`, it won't produce the filtered raw FastQ read, instead, it summarises the filtered result as UMI matrix. 

* Not to mention, it is extremely computionaly efficient. 

```
> ./fastF bam2db --help

Filter bam file with desired cell proportion and read depth with sqlite3, then summarise it UMI matrix.

    -h, --help            
    -b, --bam=<str>       path to bam file
    -f, --feature=<str>   path to feature list file
    -a, --barcode=<str>   path to barcode list file
    -d, --dbname=<str>    name of database
    -c, --cell=<flt>      rate of cell barcode
    -r, --depth=<flt>     rate of depth
    -o, --out=<str>       path to output directory
    -s, --seed=<int>      seed for random number generator
```

### `extract`: Extract the tag of bam file and write it to tag_summary.csv.

* This command just can help extract any tag from the [bam](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam) file. 

```
> ./fastF extract --help

Extract the tag of bam file.

    -h, --help        
    -b, --bam=<str>   path to bam file
    -t, --tag=<str>   tag of bam file
    -T, --type=<int>  type of tag, 0: string, 1: integer
```


