#ifndef COUNT_H
#define COUNT_H

#include "extract.h"

// get both cell barcode from R1
node *cell_counts (gzFile R1_file, size_t len_cellbarcode, size_t len_umi);

// A data structure for UMI and its collection of quality scores in different copy

typedef struct UMI_node
{
    char *UMI;
    size_t count;
    node *xf;
    struct UMI_node *left;
    struct UMI_node *right;
} UMI_node;


// A data structure for storing gene and its corresponding UMI
typedef struct gene_UMI_node
{
    char *geneID;
    UMI_node *UMI;
    size_t count;
    struct gene_UMI_node *left;
    struct gene_UMI_node *right;
} gene_UMI_node;

// A data structure for storing UMI matrix
typedef struct cell_gene_node
{
    char *CB;
    node *CR;
    gene_UMI_node *gene_UMI;
    struct cell_gene_node *left;
    struct cell_gene_node *right;
} cell_gene_node;

// new UMI_node
UMI_node *new_UMI_node(char *UMI, char *xf);

// insert UMI and geneID to UMI_node
UMI_node *insert_UMI_node(UMI_node *root, char *UMI, char *xf);

// free UMI_node
void free_UMI_node(UMI_node *root);

// print UMI_node
void print_UMI_node(UMI_node *root, gzFile fp);

// new gene_UMI_node
gene_UMI_node *new_gene_UMI_node(char *geneID, char *UMI, char *xf);

// insert a UMI to gene_UMI_node
gene_UMI_node *insert_gene_UMI_node(gene_UMI_node *root, char *geneID, char *UMI, char *xf);

// void to free the gene_UMI_node tree
void free_gene_UMI_node(gene_UMI_node *root);

// void to write the gene_UMI_node tree to a file stream;
void print_gene_UMI_node(gene_UMI_node *root, gzFile fp);

// new cell_gene_node
cell_gene_node *new_cell_gene_node(char *CB, char *CR, char *UMI, char *geneID, char *xf);

// insert a UMI to cell_gene_node
cell_gene_node *insert_cell_gene_node(cell_gene_node *root, char *CB, char *CR, char *geneID, char *UMI, char *xf);

// void to free the cell_gene_node tree
void free_cell_gene_node(cell_gene_node *root);

// count cell_gene_node
size_t count_cell_gene_node(cell_gene_node *root);

// void to write the cell_gene_node tree to a file stream;
void print_cell_gene_node(cell_gene_node *root, gzFile fp);

// count UMI nodes
size_t count_UMI_node(UMI_node *root);

// count total UMI nodes
size_t count_total_UMI_node(cell_gene_node *root);

// struct to hold the ensemble ID and its rownumber from the ensemble file

typedef struct geneID_name_node
{
    char *ensembleID;
    char *gene_name;
    struct geneID_name_node *left;
    struct geneID_name_node *right;
} geneID_name_node;

// new geneID_name_node
geneID_name_node *new_geneID_name_node(char *ensembleID, char *gene_name);

// insert ensembleID and rownumber to geneID_name_node
geneID_name_node *insert_geneID_name_node(geneID_name_node *root, char *ensembleID, char *gene_name);

// count geneID_name_node
size_t count_geneID_name_node(geneID_name_node *root);

// assume node is sorted and ensembleID is in the node
void inorder_index_geneID_name_node(geneID_name_node *root, char *ensembleID, int *index);

// search geneID_name_node
int search_geneID_for_rownumber(geneID_name_node *root, char *ensembleID);

// free geneID_name_node
void free_geneID_name_node(geneID_name_node *root);

// struct to return both cell_gene_node and geneID_name_node
typedef struct cell_gene_node_geneID_name_node
{
    cell_gene_node *cell_gene_root;
    geneID_name_node *geneID_name_root;
} cell_gene_node_geneID_name_node;

// read bam file and extract CB and gene and store them in cell_gene_node
cell_gene_node_geneID_name_node *sample_bam_UMI(
    char *bam_file, 
    char *CB_list,
    unsigned int all_cell,
    double rate_reads,
    unsigned int seed
);

// write gene_UMI_node to file stream
void write_gene_UMI_node_output(
    size_t cell_counts,
    gene_UMI_node *gene_UMI_root, 
    geneID_name_node *geneID_name_root,
    gzFile fp
);

// write cell_UMI_node to file stream
void write_cell_UMI_node_output(
    size_t *cell_counts,
    cell_gene_node *cell_gene_root,
    geneID_name_node *geneID_name_root,
    gzFile fp_barcodes,
    gzFile fp_matrix
);

// write geneID_name_node to file stream
void write_geneID_name_node_output(
    geneID_name_node *geneID_name_root,
    gzFile fp
);

#endif