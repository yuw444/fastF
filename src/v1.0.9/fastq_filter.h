// Author : Yu Wang
// Date : 01/06/2023
/*

This program is used to filter the file files based on the whitelist.
It is designed to be run with multiple threads.

*/

#ifndef QUEUE_H
#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <string.h>
#include <stdbool.h>
#include <pthread.h>

#define MAX_QUEUE_SIZE 100
#define MAX_LINE_LENGTH 256
#define LEN_CELLBARCODE 16
#define LEN_WHITELIST 2 ^ 15


int flag = 0; // flag to indicate the end of the program

pthread_mutex_t queue_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t check_lock = PTHREAD_MUTEX_INITIALIZER;

pthread_cond_t not_empty = PTHREAD_COND_INITIALIZER;
pthread_cond_t not_full = PTHREAD_COND_INITIALIZER;

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

void free_fastq(fastq *block)
{
    free(block->id);
    free(block->seq);
    free(block->qual);
    free(block);
}
void free_comb_fastq(comb_fastq *comb)
{
    free_fastq(comb->I1);
    free_fastq(comb->R1);
    free_fastq(comb->R2);

    free(comb);
}

typedef struct queue
{
    comb_fastq *data[MAX_QUEUE_SIZE];
    int head;
    int tail;
    int size;
} queue;

queue *init_queue()
{
    queue *q = (queue *)malloc(sizeof(struct queue));
    q->head = 0;
    q->tail = 0;
    q->size = 0;
    return q;
}


bool is_empty(struct queue *q)
{
    return q->size == 0;
}

bool is_full(struct queue *q)
{
    return q->size == MAX_QUEUE_SIZE;
}

void enqueue(struct queue *q, comb_fastq *value)
{
    pthread_mutex_lock(&queue_lock);
    while (is_full(q))
    {
        printf("queue is full, waiting deque!\n");
        pthread_cond_wait(&not_full, &queue_lock);
    }
    value->random_number = (double)rand() / (double)RAND_MAX;
    q->data[q->tail] = value;
    q->tail = (q->tail + 1) % MAX_QUEUE_SIZE;
    q->size++;
    pthread_cond_signal(&not_empty);
    pthread_mutex_unlock(&queue_lock);
}

comb_fastq *dequeue(struct queue *q)
{
    pthread_mutex_lock(&queue_lock);
    while (is_empty(q))
    {
        printf("queue is empty, waiting enque!\n");
        pthread_cond_wait(&not_empty, &queue_lock);
    }
    comb_fastq *value = q->data[q->head];
    q->head = (q->head + 1) % MAX_QUEUE_SIZE;
    q->size--;
    pthread_cond_signal(&not_full);
    pthread_mutex_unlock(&queue_lock);

    return value;
}

fastq *get_fastq(gzFile file)
{
    fastq *out = (fastq *)malloc(sizeof(fastq));
    out->id = (char *)malloc(MAX_LINE_LENGTH);
    out->seq = (char *)malloc(MAX_LINE_LENGTH);
    out->qual = (char *)malloc(MAX_LINE_LENGTH);

    char temp[MAX_LINE_LENGTH];
    if (gzgets(file, out->id, MAX_LINE_LENGTH) != NULL)
    {
        gzgets(file, out->seq, MAX_LINE_LENGTH);
        gzgets(file, temp, MAX_LINE_LENGTH);
        gzgets(file, out->qual, MAX_LINE_LENGTH);
    }
    else
    {
        flag = 1;
        free_fastq(out);
        return NULL;
    }

    return out;
}

comb_fastq *get_comb_fastq(gzFile fastq[3])
{
    comb_fastq *out = (comb_fastq *)malloc(sizeof(comb_fastq));
    out->I1 = get_fastq(fastq[0]);
    out->R1 = get_fastq(fastq[1]);
    out->R2 = get_fastq(fastq[2]);

    if (out->I1->id == NULL || out->R1->id == NULL || out->R2->id == NULL)
    {
        free_comb_fastq(out);
        return NULL;
    }
    return out;
}

// struct node of binary searching tree
typedef struct node
{
    char *data;
    struct node *left;
    struct node *right;
} node;

// sorting function for whitelist

int compare(const void *a, const void *b)
{
    return strcmp(*(const char **)a, *(const char **)b);
}

// construct a binary searching tree for the whitelist
node *construct_tree(char **sorted_string, int start, int end)
{
    if (start > end)
    {
        return NULL;
    }

    int mid = (start + end) / 2;
    node *root = (node *)malloc(sizeof(node));
    root->data = sorted_string[mid];
    root->left = construct_tree(sorted_string, start, mid - 1);
    root->right = construct_tree(sorted_string, mid + 1, end);

    return root;
}

// print the binary searching tree of whitelist
void print_tree(node *root)
{
    if (root == NULL)
    {
        return;
    }
    printf("%s", root->data);
    print_tree(root->left);
    print_tree(root->right);
}

// free binary tree node

void free_tree_node(node *root)
{
    if (root == NULL)
    {
        return;
    }
    free(root->data);
    free_tree_node(root->left);
    free_tree_node(root->right);
    free(root);
}

// get row number of text file

int get_row(char *file_name)
{
    int i = 0;
    char buffer[MAX_LINE_LENGTH];
    FILE *stream = fopen(file_name, "r");
    if (stream == NULL)
    {
        printf("Error opening %s!\n", file_name);
        exit(1);
    }

    while (fgets(buffer, 1024, stream) != NULL)
    {
        i++;
        // printf("%s", buffer);
    }
    fclose(stream);

    return i;
}

// read whitelist from file

char **read_txt(char *file_name, size_t nrows)
{
    FILE *stream = fopen(file_name, "r");

    if (stream == NULL)
    {
        printf("Error opening %s!\n", file_name);
        exit(1);
    }

    char **whitelist = (char **)malloc(nrows * sizeof(char *));
    int i = 0;

    for (int i = 0; i < nrows; i++)
    {
        whitelist[i] = (char *)malloc((LEN_CELLBARCODE + 2) * sizeof(char));
        fgets(whitelist[i], LEN_CELLBARCODE + 2, stream);
        // printf("%s %ld", whitelist[i], strlen(whitelist[i]));
    }

    fclose(stream);

    return whitelist;
}

bool in(node *root, char *element)
{
    if (root == NULL)
    {
        return 0;
    }

    if (strncmp(root->data, element, LEN_CELLBARCODE) == 0)
    {
        return 1;
    }
    else if (strncmp(root->data, element, LEN_CELLBARCODE) > 0)
    {
        return in(root->left, element);
    }
    else
    {
        return in(root->right, element);
    }
}

// get substring of a string

char *substring(char *string, int position, int length)
{
    char *pointer = malloc(length + 1);

    if (pointer == NULL)
    {
        printf("Unable to allocate memory.\n");
        exit(1);
    }

    strncpy(pointer, string + position, length);

    *(pointer + length) = '\0';

    return pointer;
}

// combine struct fastq to one whole string

char *combine_string(fastq *block)
{
    char *buf = (char *)malloc(MAX_LINE_LENGTH * sizeof(char));
    snprintf(buf, 512, "%s%s+\n%s", block->id, block->seq, block->qual);
    return buf;
}

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