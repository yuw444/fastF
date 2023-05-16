#include "../src/fhash/hashtable.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <htslib/sam.h>
#include <zlib.h>

#define MAX_LINE_LENGTH 1024

// struct node of binary searching tree
struct _node
{
    char *key;
    struct _node *left;
    struct _node *right;
};

typedef struct _node node;

// create a new node
node *new_node(char *key)
{
    node *newnode = (node *)malloc(sizeof(node));
    newnode->key = strdup(key);
    newnode->left = NULL;
    newnode->right = NULL;
    return newnode;
}

// insert a node to binary searching tree
node *insert_tree(node *root, char *key)
{
    if (root == NULL)
    {
        return new_node(key);
    }
    
    int cmp = strcmp(key, root->key);

    if (cmp < 0) {
        root->left = insert_tree(root->left, key);
    }
    
    if (cmp > 0) {
        root->right = insert_tree(root->right, key);
    }

    return root;
}

// free the binary searching tree
void free_tree(node *root)
{
    if (root == NULL)
    {
        return;
    }
    free(root->key);
    free_tree(root->left);
    free_tree(root->right);
    free(root);
}

// print the binary searching tree
void print_tree(node *root, FILE *stream)
{
    if (root == NULL)
    {
        return;
    }
    print_tree(root->left, stream);
    fprintf(stream, "%s\n", root->key);
    print_tree(root->right, stream);
}

size_t cell_index = 0;
// convert bst to hash table 
void bst_to_hash(node *root, hash_table **ht)
{
    if (root == NULL)
    {
        return;
    }
    bst_to_hash(root->left, ht);
    size_t *index = calloc(1, sizeof(size_t));
    *index = ++cell_index;
    hash_table_insert(*ht, root->key, index);
    bst_to_hash(root->right, ht);
}

// hash function for cell barcodes
uint64_t hash(const char *str, size_t len)
{
    uint64_t hash = 5381;
    
    for (int i=0; i < len; i++)
        hash = ((hash << 5) + hash) + str[i]; /* hash * 33 + c */

    return hash;
}


// read bam file and output cell gene UMI for each
void *write_cell_gene_UMI(
    char *bam_file,
    hash_table *ht_cell,
    hash_table *ht_gene,
    char *output_file)
{


    // FILE *output = fopen(output_file, "w");
    gzFile output = gzopen(output_file, "w");

    if (output == NULL)
    {
        printf("ERROR: Cannot open output file %s\n", output_file);
        exit(1);
    }

    // open bam file
    samFile *bam_reader = hts_open(bam_file, "r");

    if (bam_reader == NULL)
    {
        printf("ERROR: Cannot open bam file %s\n", bam_file);
        exit(1);
    }

    // read bam header
    bam_hdr_t *bam_header = sam_hdr_read(bam_reader);

    // initialize bam record to stroe each read
    bam1_t *bam_record = bam_init1();

    // initialize cell barcode tree

    size_t total_read_counts = 0;
    size_t valid_read_counts = 0;

    printf("Start processing reads\n");

    // read bam file
    while (sam_read1(bam_reader, bam_header, bam_record) >= 0)
    {
        total_read_counts++;

        // printf("total_read_counts: %lu\n", total_read_counts);

        if (total_read_counts % 1000000 == 0)
        {
            // time stamp
            time_t t = time(NULL);

            // print time stamp
            struct tm *tm = localtime(&t);
            char s[64];
            strftime(s, sizeof(s), "%c", tm);
            printf("%s: ", s);

            printf("Processed %lu reads \n", total_read_counts);
        }

        uint8_t *cb = bam_aux_get(bam_record, "CB");

        // as some read may not be aligned to any gene, we need to check if TX tag exists
        if (cb != NULL)
        {
            char *cell_barcode = bam_aux2Z(cb);
            void *temp1 = hash_table_lookup(ht_cell, cell_barcode);
            size_t cell_index;
            if(temp1 == NULL)
            {
                cell_index = 0;
            }
            else
            {
                cell_index = *(size_t *)temp1;
                valid_read_counts++;
            }
            // printf("%s\n", cell_barcode);

            uint8_t *xf = bam_aux_get(bam_record, "xf");
            int UMI_quality = bam_aux2i(xf);

            if (UMI_quality == 25)
            {
                uint8_t *gx = bam_aux_get(bam_record, "GX");
                uint8_t *ub = bam_aux_get(bam_record, "UB");
                char *gene_ID = bam_aux2Z(gx);
                void *temp2 = hash_table_lookup(ht_gene, gene_ID);
                size_t gene_index;
                if(temp2 == NULL)
                {
                    gene_index = 0;
                }
                else
                {
                    gene_index = *(size_t *)temp2;
                }
                char *UMI = bam_aux2Z(ub);
                
                gzprintf(output, "%zu,%zu,%s\n", cell_index, gene_index, UMI);

            }
        }
    }

    printf("Processed all %lu reads\n", total_read_counts);
    printf("Valid reads: %lu\n", valid_read_counts);

    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    hts_close(bam_reader);

    // fclose(output);
    gzclose(output);

}


int main(int argc, char *argv[])
{

    clock_t start = clock();

    char *file_cell = argv[1];
    char *file_gene = argv[2];
    char buffer[MAX_LINE_LENGTH];

    const int tablesize_cell = (1<<15);

    FILE *fp_cell = fopen(file_cell, "r");
    if(fp_cell == NULL)
    {
        printf("Error opening file %s\n", file_cell);
        return EXIT_FAILURE;
    }

    // read cell barcodes into bst
    printf("Start reading cell barcodes\n");
    node *root_cell = NULL;

    while(!feof(fp_cell) && fgets(buffer, MAX_LINE_LENGTH, fp_cell) != NULL)
    {
        buffer[strcspn(buffer, "\n\r\t")] = '\0';
        root_cell = insert_tree(root_cell, buffer);
    }
    
    fclose(fp_cell);

    printf("Finish reading cell barcodes\n");

    // print_tree(root_cell, stdout);
    // free_tree(root_cell);

    // save cell barcodes to hash table with object being the inorder index in bst
    // starting from 1
    
    hash_table *ht_cell = hash_table_create(tablesize_cell, hash, NULL);
    
    // convert bst to hash table
    printf("Start converting cell barcodes to hash table\n");
    bst_to_hash(root_cell, &ht_cell);
    printf("Finish converting cell barcodes to hash table\n");

    // hash_table_print(ht_cell);

    // hash table gene names
    const int tablesize_gene = (1<<16);
    hash_table *ht_gene = hash_table_create(tablesize_gene, hash, NULL);

    FILE *fp_gene = fopen(file_gene, "r");
    if(fp_gene == NULL)
    {
        printf("Error opening file %s\n", file_gene);
        return EXIT_FAILURE;
    }

    size_t gene_index = 0;
    printf("Start reading gene names\n");
    while(!feof(fp_gene) && fgets(buffer, MAX_LINE_LENGTH, fp_gene) != NULL)
    {
        buffer[strcspn(buffer, "\n\r\t")] = '\0';
        size_t *index = calloc(1, sizeof(size_t));
        *index = ++gene_index;
        hash_table_insert(ht_gene, buffer, index);
    }
    fclose(fp_gene);
    printf("Finish reading gene names\n");
    

    char *bam_file = argv[3];
    char *output_file = argv[4];

    write_cell_gene_UMI(bam_file, ht_cell, ht_gene, output_file);

    // hash_table_print(ht_cell);
    // hash_table_print(ht_gene);
    hash_table_destroy(ht_gene);
    hash_table_destroy(ht_cell);
    free_tree(root_cell);

    clock_t end = clock();

    printf("Total time: %f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);

    return 0;


}

// create a hash table to store wanted gene names