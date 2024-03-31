## fastF

* `fastF` aims to handle single-cell RNA sequencing result, FastQ file. In particular, it consists 5 subcommands. 
    * `freq`:    Find all the cell barcode whitelist and their frequencies.
    * `filter`:  Filter fastq file using cell barcode whitelist and read depth.
    * `crb`:     Extract CR and CB tags from bam file and summarize them with frequencies to a tsv file.
    * `bam2db`:  Filter bam file with desired cell proportion and read depth, then summarise it UMI matrix.
    * `extract`: Extract the tag of bam file and write it to tag_summary.csv.