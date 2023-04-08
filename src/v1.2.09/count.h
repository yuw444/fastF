#ifndef COUNT_H
#define COUNT_H

#include "extract.h"

// get both cell barcode from R1
node *cell_counts (gzFile R1_file, size_t len_cellbarcode, size_t len_umi);

// A data structure for storing UMI and its corresponding gene
typedef struct UMI_node
{
    char *UMI;
    char *geneID;
    node *xf; // a score of the UMI, as the UMI is not unique in bam file
    struct UMI_node *left;
    struct UMI_node *right;
} UMI_node;

// new UMI_node
UMI_node *new_UMI_node(char *UMI, char *geneID, char *xf);

// insert UMI and geneID to UMI_node
UMI_node *insert_UMI_node(UMI_node *root, char *UMI, char *geneID, char *xf);

// free UMI_node
void free_UMI_node(UMI_node *root);

// A data structure for storing UMI matrix
typedef struct cell_UMI_node
{
    char *CB;
    UMI_node *UMI;
    struct cell_UMI_node *left;
    struct cell_UMI_node *right;
} cell_UMI_node;

// new cell_UMI_node
cell_UMI_node *new_cell_UMI_node(char *CB, char *UMI, char *geneID, char *xf);

// insert a UMI to cell_UMI_node
cell_UMI_node *insert_cell_UMI_node(cell_UMI_node *root, char *CB, char *UMI, char *geneID, char *xf);

// void to free the cell_UMI_node tree
void free_cell_UMI_node(cell_UMI_node *root);

// void to write the cell_UMI_node tree to a file stream;
void print_cell_UMI_node(cell_UMI_node *root, FILE *fp);

// read bam file and extract CB and gene and store them in cell_UMI_node
cell_UMI_node *sample_bam_UMI(
    char *bam_file, 
    char *CB_list,
    char *gene_feature_list, 
    double rate_reads
);

#endif