
### The structure of cell_umi matrix

* `cell_gene_node` is a tree struct for both cell and its genes

    * `char *CB`: the corrected cell barcode
    * `node *CR`: the binary search tree
        
        * `char *data`: the original readed cell barcode
        * `size_t count`: the corresponding counts

    * `gene_UMI_node`: the tree struct for both gene and its corresponding UMI collection

        * `char *geneID`: the gene ID
        * `UMI_node`: the UMI collection correspond to this gene ID and UMI score collection
            * `char *UMI`: UMI
            * `node *xf`: its quality score collection 

                * `char *data`: the original readed cell barcode
                * `size_t count`: the corresponding counts


### Idea of sample bam

* Intialial run to gather the require info for downsampling

    * `node whitelist` from the reference seurat;
    * total reads `size_t total_read_counts;`
    * reads that come from cells, against the total reads count, calculate the ratio as the `amplifier ratio` later to amplifier the downsample reads `size_t valid_read_counts;` 


* Downsample
    1. Make sure the current read has a `CB` tag;
    2. Check if `CB` in `whitelist`; ratio `N` takes effect from `whitelist`;
    3. Addition condition, see if `RNG < R`; 
    4. If so, `valid_read_counts++`;
    5. Get `xf` tag;
    6. Check if `xf == 25`;
    7. If so, get `CR`, `UB`, `GX`, `GN`;
    8. `insert_cell_gene_node()`
    9. write `log.out` about `valid_read_counts++`, with the `amplifier ratio`, we could get the total raw read are need for the current alignment approximately.

* stdout
    1. Read `features.tsv.gz`, save the target `ensemble ID`
    2. process `cell_gene_node()` by cells
    3. Output both `barcode.tsv.gz` and `matrix.mtx.gz`

* format of `matrix.mtx.gz`
    1. %% is comment indicator
    2. no need to reproduce the comment content in our output
    3. But the cell number, feature counts, total feature counts are important to produce
