#include "count.h"

node *cell_counts (gzFile R1_file, size_t len_cellbarcode, size_t len_umi)
{
    node *root = NULL;

    while (1)
    {
        fastq *R1_block = get_fastq(R1_file);
        if (R1_block == NULL)
        {
            break;
        }
        char *cell_barcode_UMI = substring(R1_block->seq, 0, len_cellbarcode + len_umi);
        root = insert_tree(root, cell_barcode_UMI);
        free(cell_barcode_UMI);
        free_fastq(R1_block);
    }

    return root;
}

UMI_node *insert_UMI_node(UMI_node *root, char *CB, char *gene)
{
    if (root == NULL)
    {
        root = (UMI_node *)malloc(sizeof(UMI_node));
        root->CB = strdup(CB);
        root->gene = new_node(gene);
        root->left = NULL;
        root->right = NULL;
        return root;
    }
    else
    {
        int cmp = strcmp(CB, root->CB);
        if (cmp < 0)
        {
            root->left = insert_UMI_node(root->left, CB, gene);
        }
        else if (cmp > 0)
        {
            root->right = insert_UMI_node(root->right, CB, gene);
        }
        else
        {
            root->gene = insert_tree(root->gene, gene);
        }
    }
    return root;
}

void free_UMI_node(UMI_node *root)
{
    if (root == NULL)
    {
        return;
    }
    free(root->CB);
    free_tree_node(root->gene);

    free_UMI_node(root->left);
    free_UMI_node(root->right);
    free(root);
}

UMI_node *sample_bam_UMI(
    char *bam_file, 
    char *CB_list,
    char *gene_feature_list, 
    double rate_reads)
{
    //open CB list file
    FILE *CB_list_file = fopen(CB_list, "r");
    if (CB_list_file == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open CB list file %s\n", CB_list);
        exit(1);
    }

    // initialize CB list tree
    node *CB_list_tree = NULL;

    // read CB list file
    char CB[MAX_LINE_LENGTH];
    while (fgets(CB, MAX_LINE_LENGTH, CB_list_file) != NULL)
    {
        CB_list_tree = insert_tree(CB_list_tree, CB);
    }
    // close CB list file
    fclose(CB_list_file);

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
    UMI_node *root = NULL;

    // counter for reads
    long unsigned int read_count = 0;

    // read bam file
    while (sam_read1(bam_reader, bam_header, bam_record) >= 0){

        // extract TX tag from bam record
        uint8_t *tx = bam_aux_get(bam_record, "TX");

        // as some read may not be aligned to any gene, we need to check if TX tag exists
        if (tx != NULL)
        {
            // extract gene ID from TX tag
            char *gene = bam_aux2Z(tx);
            // string split gene ID by ","
            char *gene_ID = strtok(gene, ",");
            printf("%s\n", gene_ID);
        }
    }
    
    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    hts_close(bam_reader);

    free_tree_node(CB_list_tree);
    return root;
}
