// Author : Yu Wang
// Date : 01/06/2023
/*

This program is used to filter the file files based on the whitelist.
It is designed to be run with multiple threads.

*/

#ifndef FASTQ_FILTER_H
#define FASTQ_FILTER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#define MAX_LINE_LENGTH 128
#define LEN_CELLBARCODE 16
#define LEN_WHITELIST 2 ^ 15

// struct node of binary searching tree
typedef struct node
{
    char *data;
    struct node *left;
    struct node *right;
} node;

typedef struct fastq
{
    char *id;
    char *seq;
    char *qual;
} fastq;

typedef struct comb_fastq
{
    fastq *I1;
    fastq *R1;
    fastq *R2;
    double random_number;
} comb_fastq;

void free_fastq(fastq *block);

void free_comb_fastq(comb_fastq *comb);

fastq *get_fastq(gzFile file);

int get_comb_fastq(gzFile fastq[3], comb_fastq **block);

// sorting function for whitelist

int compare(const void *a, const void *b);

// construct a binary searching tree for the whitelist
node *construct_tree(char **sorted_string, int start, int end);

// print the binary searching tree of whitelist
void print_tree(node *root);

// free binary tree node

void free_tree_node(node *root);

// get row number of text file

int get_row(char *file_name);

// read whitelist from file

char **read_txt(char *file_name, size_t nrows);

bool in(node *root, char *element);

// get substring of a string

char *substring(char *string, int position, int length);

// combine struct fastq to one whole string

char *combine_string(fastq *block);

void fastF(gzFile file_in[3],
           gzFile file_out[3],
           node *tree_whitelist,
           unsigned int seed,
           float rate,
           bool all_cell);

// int main() {
//     struct queue *q = init_queue();
//     enqueue(q, "hello");
//     enqueue(q, "world");
//     printf("%s\n", dequeue(q));
//     printf("%s\n", dequeue(q));

//     for(int i = 0; i < 200; i++) {
//         enqueue(q, "hello");
//         printf("%s\n", dequeue(q));
//     }

//     return 0;
// }

#endif