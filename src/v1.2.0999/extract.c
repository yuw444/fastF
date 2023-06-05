
#include "extract.h"

CB_node *insert_CB_node(CB_node *root, char *CB, char *CR)
{
    if (root == NULL)
    {
        root = malloc(sizeof(CB_node));
        root->CB = strdup(CB);
        root->CR = new_node(CR);
        root->left = NULL;
        root->right = NULL;
        return root;
    }
    int cmp = strcmp(CB, root->CB);
    if (cmp < 0)
    {
        root->left = insert_CB_node(root->left, CB, CR);
    }
    else if (cmp > 0)
    {
        root->right = insert_CB_node(root->right, CB, CR);
    }
    else
    {
        // CB is equal to (*root)->CB
        root->CR = insert_tree(root->CR, CR);
    }

    return root;
}

void free_CB_node(CB_node *root)
{
    if (root == NULL)
    {
        return;
    }

    free(root->CB);
    free_tree_node(root->CR);

    free_CB_node(root->left);
    free_CB_node(root->right);
    free(root);
}

void print_CB_node(CB_node *root, gzFile fp)
{
    if (root == NULL)
    {
        return;
    }

    // print one node of CB tree into one row of csv file
    // fprintf(fp, "%s;", root->CB);
    print_tree_same_row(root->CR, fp);
    gzprintf(fp, "\n");

    print_CB_node(root->left, fp);
    print_CB_node(root->right, fp);
}

CB_node *read_bam(char *bam_file)
{
    // open bam file
    samFile *bam_reader = hts_open(bam_file, "r");

    if (bam_reader == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open bam file %s\n", bam_file);
        exit(1);
    }

    // read bam header
    bam_hdr_t *bam_header = sam_hdr_read(bam_reader);

    // initialize bam record to stroe each read
    bam1_t *bam_record = bam_init1();

    // initialize cell barcode tree
    CB_node *root = NULL;

    // counter for reads
    long unsigned int read_count = 0;

    // read bam file
    while (sam_read1(bam_reader, bam_header, bam_record) >= 0)
    {
        // extract the CB and CR tag from bam record
        uint8_t *cb = bam_aux_get(bam_record, "CB");
        uint8_t *cr = bam_aux_get(bam_record, "CR");

        // as some reads may not be able to assign CB
        if (cb != NULL)
        {
            // printf("CR: %s, CB: %s \n", bam_aux2Z(cr), bam_aux2Z(cb));
            // insert CR and CB to the CB tree
            char *CB = (char *)malloc(18 + 1);
            char *CR = (char *)malloc(16 + 1);
            strcpy(CB, bam_aux2Z(cb));
            strcpy(CR, bam_aux2Z(cr));
            root = insert_CB_node(root, CB, CR);
            free(CB);
            free(CR);
        }

        // print progress
        read_count++;
        if (read_count % 10000000 == 0)
        {
            // time stamp
            time_t t = time(NULL);

            // print time stamp
            struct tm *tm = localtime(&t);
            char s[64];
            strftime(s, sizeof(s), "%c", tm);
            printf("%s: ", s);

            printf("Processed %lu reads \n", read_count);
        }
    }

    //
    printf("Processed all %lu reads\n", read_count);
    // free bam record
    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    sam_close(bam_reader);

    return root;
}

void extract_bam(char *bam_file, const char *tag, int type)
{
    // open bam file
    samFile *bam_reader = hts_open(bam_file, "r");

    if (bam_reader == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open bam file %s\n", bam_file);
        exit(1);
    }

    // read bam header
    bam_hdr_t *bam_header = sam_hdr_read(bam_reader);

    // initialize bam record to stroe each read
    bam1_t *bam_record = bam_init1();

    // initialize cell barcode tree
    node *root = NULL;

    // counter for reads
    long total_count = 0;
    long valid_count = 0;

    // read bam file
    while (sam_read1(bam_reader, bam_header, bam_record) >= 0)
    {
        total_count++;
        // print progress
        total_count++;
        if (total_count % 10000000 == 0)
        {
            // time stamp
            time_t t = time(NULL);

            // print time stamp
            struct tm *tm = localtime(&t);
            char s[64];
            strftime(s, sizeof(s), "%c", tm);
            printf("%s: ", s);

            printf("Processed %lu reads \n", total_count);
        }

        // extract the tag

        uint8_t *tag_ptr = bam_aux_get(bam_record, tag);
        char *tag_str = NULL;
        if (type){
            tag_str = calloc(1, sizeof(int));
        }

        if (tag_ptr != NULL)
        {
            valid_count++;
            switch (type)
            {
            case 0:
                root = insert_tree(root, bam_aux2Z(tag_ptr));
                break;
            case 1:
                sprintf(tag_str, "%d", bam_aux2i(tag_ptr));
                root = insert_tree(root, tag_str);
                free(tag_str);
                break;
            }
        }
    }

    FILE *fp = fopen("tag_summary.csv", "w");
    print_tree(root, fp);
    fclose(fp);

    free_tree_node(root);

    printf("Processed all %lu reads\n", total_count);
    printf("Valid reads: %lu\n", valid_count);
    // free bam record
    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    sam_close(bam_reader);
}
