#ifndef EXTRACT_H
#define EXTRACT_H

#include <htslib/sam.h>
#include "filter.h"

// create a struct tree node to store both CR and CB tags from bam file
typedef struct CB_node {
    node *CR; // node to store CR tags
    char *CB; // corrected barcode
    struct CB_node *left;
    struct CB_node *right;
} CB_node;

// insert a CB and CR tag into the CB_node tree
CB_node *insert_CB_node(CB_node *root, char *CB, char *CR);

// void to free the CB_node tree
void free_CB_node(CB_node *root);

// void to write the CB_node tree to a file stream;
void print_CB_node(CB_node *root, gzFile fp);

// read bam file and extract CB and CR tags and store them in CB_node
CB_node *read_bam(char *bam_file);

void extract_bam(char *bam_file, const char *tag, int type);

#endif