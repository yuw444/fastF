#ifndef COUNT_H
#define COUNT_H
#include "extract.h"

extern long total_n_features;
extern long total_cells;
extern long total_genes;

// get both cell barcode from R1
node *cell_counts (gzFile R1_file, size_t len_elements, size_t len_umi);

// A data structure for UMI and its collection of quality scores in different copy
typedef struct UMI_node
{
    char *UMI;
    struct UMI_node *left;
    struct UMI_node *right;
} UMI_node;


/* A data structure for storing gene and its corresponding UMI
@nthGene: the row number of the gene in the ensemble file
@UMI: the unique UMIs of the gene
@len_UMIs: the number of unique UMIs of the gene
@left: the left child of the node
@right: the right child of the node
*/
typedef struct gene_UMI_node
{
    size_t nthGene;
    UMI_node *UMI;
    size_t len_UMIs;
    struct gene_UMI_node *left;
    struct gene_UMI_node *right;
} gene_UMI_node;

// A data structure for storing UMI matrix
typedef struct cell_gene_node
{
    size_t nthCell;
    gene_UMI_node *gene_UMI;
    struct cell_gene_node *left;
    struct cell_gene_node *right;
} cell_gene_node;

/* new UMI_node
@UMI: the UMI
@len_UMIs: the number of unique UMIs of the gene, namely the gene expression level
@return: the pointer to the new UMI_node and same time update the len_UMIs in the gene_UMI_node
*/ 

UMI_node *new_UMI_node(char *UMI);

/* insert UMI and geneID to UMI_node
@root: the root of the UMI_node
@UMI: the UMI
@len_UMIs: the number of unique UMIs of the gene, namely the gene expression level
@return: the pointer to the new UMI_node and same time update the len_UMIs in the gene_UMI_node
*/ 

UMI_node *insert_UMI_node(UMI_node *root, char *UMI);

// free UMI_node
void free_UMI_node(UMI_node *root);

/* new gene_UMI_node
@nthGene: the row number of the gene in the ensemble file
@UMI: the unique UMIs of the gene
*/ 
gene_UMI_node *new_gene_UMI_node(size_t nthGene, char *UMI);

// insert a UMI to gene_UMI_node
gene_UMI_node *insert_gene_UMI_node(gene_UMI_node *root, size_t nthGene, char *UMI);

// void to free the gene_UMI_node tree
void free_gene_UMI_node(gene_UMI_node *root);

// new cell_gene_node
cell_gene_node *new_cell_gene_node(size_t nthCell, size_t nthGene, char *UMI);

// insert a UMI to cell_gene_node
cell_gene_node *insert_cell_gene_node(cell_gene_node *root, size_t nthCell, size_t nthGene, char *UMI);
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
typedef struct list_saver_node
{
    char *element;
    size_t nthElement;
    struct list_saver_node *left;
    struct list_saver_node *right;
} list_saver_node;

// new list_saver_node
list_saver_node *new_list_saver_node(char *element, size_t nthElement);

// insert element and rownumber to list_saver_node
list_saver_node *insert_list_saver_node(list_saver_node *root, char *element, size_t nthElement);

// search element in list_saver_node and return nthElement
size_t search_list_saver_node(list_saver_node *root, char *element, size_t len_elements);

// free list_saver_node
void free_list_saver_node(list_saver_node *root);

// read bam file and extract CB and gene and store them in cell_gene_node
cell_gene_node *sample_bam_UMI(
    char *bam_file,
    char *feature_file,
    char *barcode_file,
    double rate_reads,
    unsigned int seed);

void write_gene_UMI_node_output(
    size_t nthCell,
    gene_UMI_node *gene_UMI_root,
    gzFile fp);

// write cell_UMI_node to file stream
void write_cell_gene_node_output(
    cell_gene_node *cell_gene_root,
    gzFile fp_matrix);


bool in_list(list_saver_node *root, char *element, int len_elements);

// print list_saver_node
void print_list_saver_node(list_saver_node *root);

#endif