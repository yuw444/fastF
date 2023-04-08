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

UMI_node *new_UMI_node(char *UMI, char *geneID, char *xf)
{
    if(UMI == NULL || geneID == NULL)
    {
        printf("UMI or geneID can not be NULL when creating a new UMI node!\n");
        exit(1);
    }
    UMI_node *root = (UMI_node *)malloc(sizeof(UMI_node));
    root->UMI = strdup(UMI);
    root->geneID = strdup(geneID);
    root->xf = new_node(xf);
    root->left = NULL;
    root->right = NULL;
    return root;
}

UMI_node *insert_UMI_node(UMI_node *root, char *UMI, char *geneID, char *xf)
{
    if (root == NULL)
    {
        root = new_UMI_node(UMI, geneID, xf);
        return root;
    }
    else
    {
        int cmp = strcmp(UMI, root->UMI);
        if (cmp < 0)
        {
            root->left = insert_UMI_node(root->left, UMI, geneID, xf);
        }
        else if (cmp > 0)
        {
            root->right = insert_UMI_node(root->right, UMI, geneID, xf);
        }
        else
        {
            root->xf = insert_tree(root->xf, xf);
        }
        return root;
    }
}

void free_UMI_node(UMI_node *root)
{
    if (root == NULL)
    {
        return;
    }
    free(root->UMI);
    free(root->geneID);

    free_UMI_node(root->left);
    free_UMI_node(root->right);
    free(root);
}

cell_UMI_node *new_cell_UMI_node(char *CB, char *UMI, char *geneID, char *xf)
{
    if(CB == NULL || UMI == NULL || geneID == NULL || xf == NULL)
    {
        printf("CB or UMI or geneID or xf can not be NULL when creating a new cell_UMI node!\n");
        exit(1);
    }
    cell_UMI_node *root = (cell_UMI_node *)malloc(sizeof(cell_UMI_node));
    root->CB = strdup(CB);
    root->UMI = new_UMI_node(UMI, geneID, xf);
    root->left = NULL;
    root->right = NULL;
    return root;
}

cell_UMI_node *insert_cell_UMI_node(cell_UMI_node *root, char *CB, char *UMI, char *geneID, char *xf)
{
    if (root == NULL)
    {
        root = new_cell_UMI_node(CB, UMI, geneID, xf);
        return root;
    }
    else
    {
        int cmp = strcmp(CB, root->CB);
        if (cmp < 0)
        {
            root->left = insert_cell_UMI_node(root->left, CB, UMI, geneID, xf);
        }
        else if (cmp > 0)
        {
            root->right = insert_cell_UMI_node(root->right, CB, UMI, geneID, xf);
        }
        else
        {
            root->UMI = insert_UMI_node(root->UMI, UMI, geneID, xf);
        }
    }
    return root;
}

void free_cell_UMI_node(cell_UMI_node *root)
{
    if (root == NULL)
    {
        return;
    }
    free(root->CB);
    free_UMI_node(root->UMI);

    free_cell_UMI_node(root->left);
    free_cell_UMI_node(root->right);
    free(root);
}

cell_UMI_node *sample_bam_UMI(
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
        // get the first 16 characters as cell barcode
        CB[16] = '\0';
        CB_list_tree = insert_tree(CB_list_tree, CB);
    }
    // close CB list file
    fclose(CB_list_file);

    // open a file to store CB list tree
    FILE *CB_list_tree_file = fopen("CB_list_tree.txt", "w");
    print_tree(CB_list_tree, CB_list_tree_file);
    fclose(CB_list_tree_file);

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
    cell_UMI_node *root = NULL;

    // counter for reads
    long unsigned int read_count = 0;

    // read bam file
    while (sam_read1(bam_reader, bam_header, bam_record) >= 0){

        // extract CB, xf, GX, ub tag from bam record
        uint8_t *cb = bam_aux_get(bam_record, "CB");
        uint8_t *xf = bam_aux_get(bam_record, "xf");
        uint8_t *gx = bam_aux_get(bam_record, "GX");
        uint8_t *ub = bam_aux_get(bam_record, "UB");

        // printf("i = %ld\n", ++read_count);

        // as some read may not be aligned to any gene, we need to check if TX tag exists
        if (cb != NULL && xf != NULL && gx != NULL && ub != NULL)
        {
            // extract gene ID from TX tag
            char *cell_barcode = bam_aux2Z(cb);
            int UMI_counting = bam_aux2i(xf);
            char *gene_ID = bam_aux2Z(gx);
            char *UMI = bam_aux2Z(ub);
            if(UMI_counting == 25)
                printf("%s; %ld; %s; %s\n", cell_barcode, UMI_counting, gene_ID, UMI);
        }
    }
    
    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    hts_close(bam_reader);

    free_tree_node(CB_list_tree);
    return root;
}
