#ifndef DATA_H
#define DATA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <htslib/sam.h>
#include <htslib/hts.h>


/* UMI node
*/

typedef struct umi_node {
    char *umi;
    struct umi_node *left;
    struct umi_node *right;
} umi_node;

/* gene UMI node
*/

typedef struct gene_umi_node {
    char *gene;
    size_t nth_gene;
    size_t expression_level;
    umi_node *umi_root;
    struct gene_umi_node *left;
    struct gene_umi_node *right;
} gene_umi_node;

/* cell gene UMI node
*/

typedef struct cell_gene_umi_node {
    char *cell;
    gene_umi_node *gene_root;
    struct cell_gene_umi_node *left;
    struct cell_gene_umi_node *right;
} cell_gene_umi_node;

/* new umi node
*/

umi_node *new_umi_node(char *umi);
umi_node *insert_umi_node(umi_node *root, char *umi);
void free_umi_node(umi_node *root);

/* new gene UMI node
*/

gene_umi_node *new_gene_umi_node(char *gene, char *umi);
gene_umi_node *insert_gene_umi_node(gene_umi_node *root, char *gene, char *umi);
void free_gene_umi_node(gene_umi_node *root);

/* new cell gene UMI node
*/

cell_gene_umi_node *new_cell_gene_umi_node(char *cell, char *gene, char *umi);
cell_gene_umi_node *insert_cell_gene_umi_node(cell_gene_umi_node *root, char *cell, char *gene, char *umi);
void free_cell_gene_umi_node(cell_gene_umi_node *root);

#endif