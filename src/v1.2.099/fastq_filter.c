// Author : Yu Wang
// Date : 01/13/2023
/*

This program is used to filter the file files based on the whitelist.
It is designed to be run with multiple threads.

*/

#include "fastq_filter.h"

int flag = 0;

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

int get_comb_fastq(gzFile fastq[3], comb_fastq **block)
{
    if(fastq[0] != Z_NULL)
    {
        (*block)->I1 = get_fastq(fastq[0]);
    }
    // R1 cannot be NULL
    (*block)->R1 = get_fastq(fastq[1]);

    if(fastq[2] != Z_NULL)
    {
        (*block)->R2 = get_fastq(fastq[2]);
    }
    (*block)->random_number = (float) rand() / RAND_MAX;

    if ((*block)->R1 == NULL)
    {
        return -1;
    }
    return 0;
}
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

node *new_node(char *data)
{
    if (data == NULL)
    {
        return NULL;
    }

    node *out = (node *)malloc(sizeof(node));
    out->data = strdup(data);
    out->count = 1;
    out->left = NULL;
    out->right = NULL;
    return out;
}

void get_node_inorder_index(node *root, char *data, size_t len_data, int *index)
{
    if (root == NULL)
    {
        return;
    }
    get_node_inorder_index(root->left, data, len_data, index);
    if (strncmp(root->data, data, len_data) <= 0)
    {
        (*index)++;
    }
    get_node_inorder_index(root->right, data, len_data, index);
}

node *insert_tree(node *root, char *data)
{
    if (root == NULL)
    {
        return new_node(data);
    }
    
    int cmp = strcmp(data, root->data);

    if (cmp < 0) {
        root->left = insert_tree(root->left, data);
    }
    else if (cmp > 0) {
        root->right = insert_tree(root->right, data);
    }
    else {  // data is equal to root->data
        root->count++;
    }
    return root;
}

// construct a binary searching tree for the whitelist
node *construct_tree(char **string_vector, size_t len)
{
    node *root = NULL;
    for (int i = 0; i < len; i++)
    {
        root = insert_tree(root, string_vector[i]);
    }
    return root;
}


// print the binary searching tree of whitelist
void print_tree(node *root, FILE *stream)
{
    if (root == NULL)
    {
        return;
    }
    fprintf(stream, "%s,%ld\n", root->data, root->count);
    print_tree(root->left, stream);
    print_tree(root->right, stream);
}

void print_tree_gz(node *root, gzFile stream)
{
    if (root == NULL)
    {
        return;
    }
    print_tree_gz(root->left, stream);
    gzprintf(stream, "%s\n", root->data);
    print_tree_gz(root->right, stream);
}

void print_tree_same_row(node *root, gzFile stream)
{
    if (root == NULL)
    {
        return;
    }
    gzprintf(stream, "%s,%ld;", root->data, root->count);
    print_tree_same_row(root->left, stream);
    print_tree_same_row(root->right, stream);
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

    while (fgets(buffer, MAX_LINE_LENGTH, stream) != NULL)
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
        whitelist[i] = (char *)malloc(100 * sizeof(char));
        fgets(whitelist[i], 100, stream);
    }

    fclose(stream);

    return whitelist;
}

bool in(node *root, char *element, int len_cellbarcode)
{
    if (root == NULL)
    {
        return 0;
    }

    if (strncmp(root->data, element, len_cellbarcode) == 0)
    {
        return 1;
    }
    else if (strncmp(root->data, element, len_cellbarcode) > 0)
    {
        return in(root->left, element, len_cellbarcode);
    }
    else
    {
        return in(root->right, element, len_cellbarcode);
    }
}

// get substring of a string

char *substring(char *string, int position, int length)
{
    char *pointer = (char *)malloc(length + 1);

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
    char *buf = (char *)malloc(MAX_LINE_LENGTH * 3 * sizeof(char));
    snprintf(buf, MAX_LINE_LENGTH * 3, "%s%s+\n%s", block->id, block->seq, block->qual);
    return buf;
}

void fastF(gzFile file_in[3],
           gzFile file_out[3],
           node *tree_whitelist,
           unsigned int len_cellbarcode,
           unsigned int seed,
           float rate,
           bool  all_cell)
{
    srand(seed);

    comb_fastq *block = (comb_fastq *)malloc(sizeof(comb_fastq));
    block->I1 = NULL;
    block->R1 = NULL;
    block->R2 = NULL;

#pragma omp parallel
    {
#pragma omp single
        {
            while (get_comb_fastq(file_in, &block) != -1)
            {
#pragma omp task firstprivate(block)
                {
                    char *cell_barcode = substring(block->R1->seq, 0, len_cellbarcode);

                    if (block->random_number < rate && (all_cell || in(tree_whitelist, cell_barcode, len_cellbarcode)))
                    {
#pragma omp critical
                        {
                            if (block->I1 != NULL)
                            {
                                char *I1 = combine_string(block->I1);
                                gzputs(file_out[0], I1);
                                free(I1);
                            }

                            // R1 can not be NULL
                            char *R1 = combine_string(block->R1);
                            gzputs(file_out[1], R1);
                            free(R1);

                            if (block->R2 != NULL)
                            {
                                char *R2 = combine_string(block->R2);
                                gzputs(file_out[2], R2);
                                free(R2);
                            }
                        }
                    }
                    free(cell_barcode);
                    if(block->I1 != NULL)
                        free_fastq(block->I1);
                        
                    // R1 can not be NULL
                    free_fastq(block->R1);

                    if (block->R2 != NULL)
                        free_fastq(block->R2);
                }
            }
        }
    }

    free(block);
}

// int main()
// {
//     double startTime = clock();

//     unsigned int seed = 926;
//     srand(seed);

//     file_in[0] = gzopen("/home/rstudio/Frag/data/I1.fastq.gz", "rb");
//     file_in[1] = gzopen("/home/rstudio/Frag/data/R1.fastq.gz", "rb");
//     file_in[2] = gzopen("/home/rstudio/Frag/data/R2.fastq.gz", "rb");

//     file_out[0] = gzopen("/home/rstudio/Frag/data/I1_openmp.fastq.gz", "wb");
//     file_out[1] = gzopen("/home/rstudio/Frag/data/R1_openmp.fastq.gz", "wb");
//     file_out[2] = gzopen("/home/rstudio/Frag/data/R2_openmp.fastq.gz", "wb");

//     rate_threshold = 0.8;

//     char *filename_whitelist = "../data/whitelist.txt";

//     printf("Reading whitelist...\n");
//     int nrow = get_row(filename_whitelist);

//     printf("nrow = %d\n", nrow);

//     char **whitelist = read_txt(filename_whitelist, nrow);

//     qsort(whitelist, nrow, sizeof(char *), compare);

//     tree_whitelist = construct_tree(whitelist, 0, nrow - 1);

//     // print_tree(tree_whitelist);

//     int t1 = in(tree_whitelist, "TTTGGTTGTGACCAAG");

//     printf("t1 = %d\n", t1);

//     fastF(file_in, file_out, tree_whitelist, seed, rate_threshold);

//     free_tree_node(tree_whitelist);
//     free(whitelist);

//     for (int i = 0; i < 3; i++)
//     {
//         gzclose(file_in[i]);
//         gzclose(file_out[i]);
//     }

//     double stopTime = clock();

//     double secsElapsed = (stopTime - startTime) / CLOCKS_PER_SEC;

//     printf("Elapsed: %f seconds\n", secsElapsed);
//     return 0;
// }

// // Elapsed: 403.236648 seconds