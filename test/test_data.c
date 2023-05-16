#include "../src/bamonly/data.h"
#include <criterion/criterion.h>

Test(data, SIZE)
{
    printf("Size of cell_gene_node: %lu\n", sizeof(cell_gene_umi_node));
}

