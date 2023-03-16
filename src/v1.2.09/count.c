#include "count.h"

node *cell_counts (gzFile R1_file, size_t len_cellbarcode, size_t len_umi)
{
    node *root = NULL;

    while (1)
    {
        fastq *R1_block = get_fastq(R1_file);
        if (R1_block == NULL)
        {
            break;
        }
        char *cell_barcode_UMI = substring(R1_block->seq, 0, len_cellbarcode + len_umi);
        root = insert_tree(root, cell_barcode_UMI);
        free(cell_barcode_UMI);
        free_fastq(R1_block);
    }

    return root;
}

