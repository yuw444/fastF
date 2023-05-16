#include "../src/fhash/hashtable.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <htslib/sam.h>
#include <zlib.h>
#include <sqlite3.h>
#define MAX_LINE_LENGTH 1024
#define BASE_BITS 2                      // the number of bits to encode a base
#define BYTE_SIZE 8                      // the number of bits in a byte
#define BYTE_MASK 0xFF                   // the mask of a byte
#define BASE_MASK ((1 << BASE_BITS) - 1) // the mask of a base

sqlite3 *db;
char *zErrMsg = 0;
int rc;
char *sql;

uint8_t encode_base(char base)
{
    switch (base)
    {
    case 'A':
        return 0b00;
    case 'C':
        return 0b01;
    case 'G':
        return 0b10;
    case 'T':
        return 0b11;
    default:
        return BYTE_MASK; // unknown base
    }
}

uint8_t *encode_UMI(const char *UMI)
{

    if (strlen(UMI) != 10)
    {
        printf("Error: UMI length is not 10\n");
        exit(1);
    }
    uint8_t *encoded_UMI = (uint8_t *)calloc(4, sizeof(uint8_t));
    if (encoded_UMI == NULL)
    {
        return NULL;
    }
    uint8_t *p = encoded_UMI;
    int bit_index = BYTE_SIZE - BASE_BITS;
    for (size_t i = 0; i < 10; i++)
    {
        uint8_t base = encode_base(UMI[i]);
        if (base == BYTE_MASK)
        {
            free(encoded_UMI);
            return NULL; // unknown base
        }
        *p |= (base & BASE_MASK) << bit_index;
        bit_index -= BASE_BITS;
        if (bit_index < 0)
        {
            bit_index = BYTE_SIZE - BASE_BITS;
            p++;
        }
    }
    return encoded_UMI;
}

// decode to char* with length n from unit8_t

char *decode_DNA(uint8_t *encoded_UMI, size_t n)
{
    char *decoded_DNA = (char *)calloc(n + 1, sizeof(char));
    if (decoded_DNA == NULL)
    {
        return NULL;
    }
    char *p = decoded_DNA;
    int bit_index = BYTE_SIZE - BASE_BITS;
    for (size_t i = 0; i < n; i++)
    {
        uint8_t base = (encoded_UMI[i / 4] >> bit_index) & BASE_MASK;
        switch (base)
        {
        case 0b00:
            *p = 'A';
            break;
        case 0b01:
            *p = 'C';
            break;
        case 0b10:
            *p = 'G';
            break;
        case 0b11:
            *p = 'T';
            break;
        default:
            free(decoded_DNA);
            return NULL; // unknown base
        }
        bit_index -= BASE_BITS;
        if (bit_index < 0)
        {
            bit_index = BYTE_SIZE - BASE_BITS;
            p++;
        }
    }
    return decoded_DNA;
}

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

    if (cmp < 0)
    {
        root->left = insert_tree(root->left, key);
    }

    if (cmp > 0)
    {
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

// inorder traversal of binary searching tree, insert cell to sql table
size_t cell_index_sql = 1;
void inorder_traversal_sql_insert(node *root, sqlite3_stmt *stmt)
{
    if (root == NULL)
    {
        return;
    }
    inorder_traversal_sql_insert(root->left, stmt);

    sqlite3_bind_int(stmt, 1, cell_index_sql++);
    sqlite3_bind_text(stmt, 2, root->key, strlen(root->key), SQLITE_STATIC);

    if(sqlite3_step(stmt) != SQLITE_DONE){
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
    }
    sqlite3_reset(stmt);

    inorder_traversal_sql_insert(root->right, stmt);
}

// hash function for cell barcodes
uint64_t hash(const char *str, size_t len)
{
    uint64_t hash = 5381;

    for (int i = 0; i < len; i++)
        hash = ((hash << 5) + hash) + str[i]; /* hash * 33 + c */

    return hash;
}

// read bam file and output cell gene UMI for each
void *write_cell_gene_UMI(
    char *bam_file,
    hash_table *ht_cell,
    hash_table *ht_gene)
{

    sql = "CREATE TABLE umi (cell_index INTEGER, gene_index INTEGER, encoded_UMI TEXT);";
    rc = sqlite3_exec(db, sql, NULL, 0, &zErrMsg);

    if (rc != SQLITE_OK)
    {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
    }

    // open bam file
    samFile *bam_reader = hts_open(bam_file, "r");
    if (bam_reader == NULL)
    {
        printf("ERROR: Cannot open bam file %s\n", bam_file);
        exit(1);
    }

    // counters
    size_t total_read_counts = 0;
    size_t valid_read_counts = 0;

    // read bam header
    bam_hdr_t *bam_header = sam_hdr_read(bam_reader);
    bam1_t *bam_record = bam_init1();

    // read bam file
    printf("Start processing reads...\n");

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);
    sql = "INSERT INTO umi VALUES (?1, ?2, ?3);";
    sqlite3_stmt *stmt;
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    while (sam_read1(bam_reader, bam_header, bam_record) >= 0)
    {
        total_read_counts++;

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
            if (temp1 != NULL)
            {
                cell_index = *(size_t *)temp1;
                uint8_t *xf = bam_aux_get(bam_record, "xf");
                int UMI_quality = bam_aux2i(xf);

                if (UMI_quality == 25)
                {
                    uint8_t *gx = bam_aux_get(bam_record, "GX");
                    uint8_t *ub = bam_aux_get(bam_record, "UB");
                    char *gene_ID = bam_aux2Z(gx);
                    void *temp2 = hash_table_lookup(ht_gene, gene_ID);
                    size_t gene_index;
                    if (temp2 != NULL)
                    {
                        valid_read_counts++;
                        gene_index = *(size_t *)temp2;
                        char *UMI = bam_aux2Z(ub);
                        char *encoded_UMI = encode_UMI(UMI);
                        sqlite3_bind_int(stmt, 1, cell_index);
                        sqlite3_bind_int(stmt, 2, gene_index);
                        sqlite3_bind_text(stmt, 3, encoded_UMI, strlen(encoded_UMI), SQLITE_STATIC);
                        free(encoded_UMI);

                        if (sqlite3_step(stmt) != SQLITE_DONE)
                        {
                            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
                        }

                        sqlite3_reset(stmt);
                    }
                }
            }
        }
    }

    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);
    sqlite3_finalize(stmt);

    printf("Processed all %lu reads\n", total_read_counts);
    printf("Valid reads: %lu\n", valid_read_counts);

    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    hts_close(bam_reader);

}

int main(int argc, char *argv[])
{

    clock_t start = clock();

    char *file_cell = argv[1];
    char *file_gene = argv[2];
    char *bam_file = argv[3];
    char *name_db = argv[4];

    char buffer[MAX_LINE_LENGTH];

    const int tablesize_cell = (1 << 15);

    FILE *fp_cell = fopen(file_cell, "r");
    if (fp_cell == NULL)
    {
        printf("Error opening file %s\n", file_cell);
        return EXIT_FAILURE;
    }

    // read cell barcodes into bst
    // for the purpose of ordering
    printf("Start reading cell barcodes\n");
    node *root_cell = NULL;

    while (!feof(fp_cell) && fgets(buffer, MAX_LINE_LENGTH, fp_cell) != NULL)
    {
        buffer[strcspn(buffer, "\n\r\t")] = '\0';
        root_cell = insert_tree(root_cell, buffer);
    }

    fclose(fp_cell);
    printf("Finish reading cell barcodes\n");

    // open database
    rc = sqlite3_open("umi.db", &db);
    if (rc)
    {
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        return (1);
    }
    else
    {
        printf("Opened database successfully\n");
    }

    /*****************create table cell*****************/

    sql = "CREATE TABLE cell(" \
          "id INTEGER PRIMARY KEY AUTOINCREMENT," \
          "barcode TEXT NOT NULL);";

    rc = sqlite3_exec(db, sql, NULL, NULL, &zErrMsg);

    if(rc != SQLITE_OK)
    {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }
    else
    {
        printf("Table cell created successfully\n");
    }

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);
    sql = "INSERT INTO cell VALUES (?1, ?2);";
    sqlite3_stmt *stmt;
    sqlite3_prepare_v2(db, sql, strlen(sql), &stmt, NULL);

    printf("Start inserting cell barcodes\n");
    inorder_traversal_sql_insert(root_cell, stmt);

    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);
    sqlite3_finalize(stmt);
    printf("Finish inserting cell barcodes\n");


    /*****************create table gene*****************/

    sql = "CREATE TABLE gene(" \
          "id INTEGER PRIMARY KEY AUTOINCREMENT," \
          "name TEXT NOT NULL);";

    sqlite3_close(db);

    return 0;

    hash_table *ht_cell = hash_table_create(tablesize_cell, hash, NULL);

    // convert bst to hash table
    printf("Start converting cell barcodes to hash table\n");
    bst_to_hash(root_cell, &ht_cell);
    printf("Finish converting cell barcodes to hash table\n");

    // hash_table_print(ht_cell);

    // hash table gene names
    const int tablesize_gene = (1 << 16);
    hash_table *ht_gene = hash_table_create(tablesize_gene, hash, NULL);

    FILE *fp_gene = fopen(file_gene, "r");
    if (fp_gene == NULL)
    {
        printf("Error opening file %s\n", file_gene);
        return EXIT_FAILURE;
    }

    size_t gene_index = 0;
    printf("Start reading gene names\n");
    while (!feof(fp_gene) && fgets(buffer, MAX_LINE_LENGTH, fp_gene) != NULL)
    {
        buffer[strcspn(buffer, "\n\r\t")] = '\0';
        size_t *index = calloc(1, sizeof(size_t));
        *index = ++gene_index;
        hash_table_insert(ht_gene, buffer, index);
    }
    fclose(fp_gene);
    printf("Finish reading gene names\n");

    char *output_file = argv[4];

    write_cell_gene_UMI(bam_file, ht_cell, ht_gene);

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