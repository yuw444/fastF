#ifndef COUNT_H
#define COUNT_H

#include "extract.h"

// get both cell barcode from R1
node *cell_counts (gzFile R1_file, size_t len_cellbarcode, size_t len_umi);

// A data structure for storing UMI matrix
typedef struct UMI_node
{
    char *CB;
    node *gene;
    struct UMI_node *left;
    struct UMI_node *right;
} UMI_node;

// insert a UMI node to UMI tree
UMI_node *insert_UMI_node(UMI_node *root, char *CB, char *gene);

// void to free the UMI_node tree
void free_UMI_node(UMI_node *root);

// void to write the UMI_node tree to a file stream;
void print_UMI_node(UMI_node *root, FILE *fp);

// read bam file and extract CB and gene and store them in UMI_node
UMI_node *sample_bam_UMI(
    char *bam_file, 
    char *CB_list,
    char *gene_feature_list, 
    double rate_reads
);

#endif