#include "count.h"

node *cell_counts(gzFile R1_file, size_t len_cellbarcode, size_t len_umi)
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

UMI_node *new_UMI_node(char *UMI, char *xf)
{
    if (UMI == NULL || xf == NULL)
    {
        return NULL;
    }
    UMI_node *root = (UMI_node *)malloc(sizeof(UMI_node));
    root->UMI = strdup(UMI);
    root->xf = new_node(xf);
    root->count = 1;
    root->left = NULL;
    root->right = NULL;
    return root;
}

UMI_node *insert_UMI_node(UMI_node *root, char *UMI, char *xf)
{
    if (root == NULL)
    {
        root = new_UMI_node(UMI, xf);
        return root;
    }
    else
    {
        int cmp = strcmp(UMI, root->UMI);
        if (cmp < 0)
        {
            root->left = insert_UMI_node(root->left, UMI, xf);
        }
        else if (cmp > 0)
        {
            root->right = insert_UMI_node(root->right, UMI, xf);
        }
        else
        {
            root->xf = insert_tree(root->xf, xf);
        }
        return root;
    }
}

size_t count_UMI_node(UMI_node *root)
{
    if (root == NULL)
    {
        return 0;
    }
    else
    {
        return 1 + count_UMI_node(root->left) + count_UMI_node(root->right);
    }
}

void free_UMI_node(UMI_node *root)
{
    if (root == NULL)
    {
        return;
    }

    if (root->UMI != NULL)
    {
        free(root->UMI);
    }

    if (root->xf != NULL)
    {
        free_tree_node(root->xf);
    }

    free_UMI_node(root->left);
    free_UMI_node(root->right);
    free(root);
}

void print_UMI_node(UMI_node *root, gzFile fp)
{
    if (root == NULL)
    {
        return;
    }

    print_UMI_node(root->left, fp);
    gzprintf(fp, " %s:", root->UMI);
    print_tree_same_row(root->xf, fp);
    print_UMI_node(root->right, fp);
}

gene_UMI_node *new_gene_UMI_node(char *geneID, char *UMI, char *xf)
{
    if (geneID == NULL || UMI == NULL || xf == NULL)
    {
        return NULL;
    }
    gene_UMI_node *root = (gene_UMI_node *)malloc(sizeof(gene_UMI_node));
    root->geneID = strdup(geneID);
    root->UMI = new_UMI_node(UMI, xf);
    root->left = NULL;
    root->right = NULL;
    return root;
}

gene_UMI_node *insert_gene_UMI_node(
    gene_UMI_node *root,
    char *geneID,
    char *UMI,
    char *xf)
{
    if (root == NULL)
    {
        root = new_gene_UMI_node(geneID, UMI, xf);
        return root;
    }
    else
    {
        int cmp = strcmp(geneID, root->geneID);
        if (cmp < 0)
        {
            root->left = insert_gene_UMI_node(root->left, geneID, UMI, xf);
        }
        else if (cmp > 0)
        {
            root->right = insert_gene_UMI_node(root->right, geneID, UMI, xf);
        }
        else
        {
            root->UMI = insert_UMI_node(root->UMI, UMI, xf);
        }
        return root;
    }
}

void free_gene_UMI_node(gene_UMI_node *root)
{
    if (root == NULL)
    {
        return;
    }

    if (root->geneID != NULL)
    {
        free(root->geneID);
    }

    if (root->UMI != NULL)
    {
        free_UMI_node(root->UMI);
    }

    free_gene_UMI_node(root->left);
    free_gene_UMI_node(root->right);

    free(root);
}

void print_gene_UMI_node(gene_UMI_node *root, gzFile fp)
{
    if (root == NULL)
    {
        return;
    }

    print_gene_UMI_node(root->left, fp);
    gzprintf(fp, "\t %s:", root->geneID);
    print_UMI_node(root->UMI, fp);
    print_gene_UMI_node(root->right, fp);
}

cell_gene_node *new_cell_gene_node(char *CB, char *CR, char *geneID, char *UMI, char *xf)
{
    if (CB == NULL || UMI == NULL || geneID == NULL || xf == NULL)
    {
        printf("CB or UMI or geneID or xf can not be NULL when creating a new cell_gene_node!\n");
        exit(1);
    }
    cell_gene_node *root = (cell_gene_node *)malloc(sizeof(cell_gene_node));
    root->CB = strdup(CB);
    root->CR = new_node(CR);
    root->gene_UMI = new_gene_UMI_node(geneID, UMI, xf);
    root->left = NULL;
    root->right = NULL;
    return root;
}

cell_gene_node *insert_cell_gene_node(cell_gene_node *root, char *CB, char *CR, char *geneID, char *UMI, char *xf)
{
    if (root == NULL)
    {
        root = new_cell_gene_node(CB, CR, geneID, UMI, xf);
        return root;
    }
    else
    {
        int cmp = strcmp(CB, root->CB);
        if (cmp < 0)
        {
            root->left = insert_cell_gene_node(root->left, CB, CR, geneID, UMI, xf);
        }
        else if (cmp > 0)
        {
            root->right = insert_cell_gene_node(root->right, CB, CR, geneID, UMI, xf);
        }
        else
        {
            root->gene_UMI = insert_gene_UMI_node(root->gene_UMI, geneID, UMI, xf);
        }
    }
    return root;
}

void free_cell_gene_node(cell_gene_node *root)
{
    if (root == NULL)
    {
        return;
    }
    if (root->CB != NULL)
    {
        free(root->CB);
    }

    if (root->CR != NULL)
    {
        free_tree_node(root->CR);
    }

    if (root->gene_UMI != NULL)
    {
        free_gene_UMI_node(root->gene_UMI);
    }

    free_cell_gene_node(root->left);
    free_cell_gene_node(root->right);

    free(root);
}

void print_cell_gene_node(cell_gene_node *root, gzFile fp)
{
    if (root == NULL)
    {
        return;
    }

    print_cell_gene_node(root->left, fp);
    print_tree_same_row(root->CR, fp);
    print_gene_UMI_node(root->gene_UMI, fp);
    gzprintf(fp, "\n");
    print_cell_gene_node(root->right, fp);
}

size_t count_cell_gene_node(cell_gene_node *root)
{
    if (root == NULL)
    {
        return 0;
    }
    return 1 + count_cell_gene_node(root->left) + count_cell_gene_node(root->right);
}

size_t count_total_UMI_node(cell_gene_node *root)
{
    if (root == NULL)
    {
        return 0;
    }
    return count_UMI_node(root->gene_UMI->UMI) + count_total_UMI_node(root->left) + count_total_UMI_node(root->right);
}

cell_gene_node_geneID_name_node *sample_bam_UMI(
    char *bam_file,
    char *CB_list,
    unsigned int all_cell,
    double rate_reads,
    unsigned int seed)
{
    srand(seed);

    // initialize CB list tree
    node *CB_list_tree = NULL;

    if (all_cell == 0)
    {
        // open CB list file
        FILE *CB_list_file = fopen(CB_list, "r");
        if (CB_list_file == NULL)
        {
            printf("ERROR: Cannot open CB list file %s\n", CB_list);
            exit(1);
        }
        
        // read CB list file
        char CB[MAX_LINE_LENGTH];
        while (fgets(CB, MAX_LINE_LENGTH, CB_list_file) != NULL)
        {
            // remove newline character
            CB[16] = '\0';
            CB_list_tree = insert_tree(CB_list_tree, CB);
        }
        // close CB list file
        fclose(CB_list_file);
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
    cell_gene_node *cell_gene_root = NULL;

    // initialize gene ID tree
    geneID_name_node *geneID_name_root = NULL;

    size_t total_read_counts = 0;

    // counter for reads
    size_t valid_read_counts = 0;

    // read bam file
    while (sam_read1(bam_reader, bam_header, bam_record) >= 0)
    {
        total_read_counts++;

        // print progress
        total_read_counts++;

        if (total_read_counts % 10000000 == 0)
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

        // extract CB, xf, GX, ub tag from bam record
        uint8_t *cb = bam_aux_get(bam_record, "CB");

        // get a random number between 0 and 1
        double random_number = (double)rand() / (double)RAND_MAX;

        // as some read may not be aligned to any gene, we need to check if TX tag exists
        if (cb != NULL)
        {
            // extract gene ID from TX tag
            char *cell_barcode = bam_aux2Z(cb);
            // printf("%s\n", cell_barcode);

            if (all_cell || in(CB_list_tree, cell_barcode, 16))
            {
                // printf("%s\n", cell_barcode);

                if (random_number < rate_reads)
                {
                    valid_read_counts++;
                    // get xf tag
                    uint8_t *xf = bam_aux_get(bam_record, "xf");
                    int UMI_quality = bam_aux2i(xf);
                    char quality_string[3];
                    sprintf(quality_string, "%d", UMI_quality);

                    if (UMI_quality == 25)
                    {
                        uint8_t *cr = bam_aux_get(bam_record, "CR");
                        uint8_t *gx = bam_aux_get(bam_record, "GX");
                        uint8_t *gn = bam_aux_get(bam_record, "GN");
                        uint8_t *ub = bam_aux_get(bam_record, "UB");
                        char *CR = bam_aux2Z(cr);
                        char *gene_ID = bam_aux2Z(gx);
                        char *gene_name = bam_aux2Z(gn);
                        char *UMI = bam_aux2Z(ub);

                        printf("%s; %ld; %s; %s\n", cell_barcode, UMI_quality, gene_ID, UMI);

                        cell_gene_root = insert_cell_gene_node(
                            cell_gene_root,
                            cell_barcode,
                            CR,
                            gene_ID,
                            UMI,
                            quality_string);

                        geneID_name_root = insert_geneID_name_node(
                            geneID_name_root,
                            gene_ID,
                            gene_name);
                    }
                }
            }
        }
    }

    printf("Processed all %lu reads\n", total_read_counts);

    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    hts_close(bam_reader);

    printf("Total reads: %zu\n", total_read_counts);
    printf("Valid reads: %zu\n", valid_read_counts);

    if (CB_list_tree != NULL)
    {
        free_tree_node(CB_list_tree);
    }

    cell_gene_node_geneID_name_node *root = (cell_gene_node_geneID_name_node *)malloc(sizeof(cell_gene_node_geneID_name_node));
    root->cell_gene_root = cell_gene_root;
    root->geneID_name_root = geneID_name_root;

    return root;
}

geneID_name_node *new_geneID_name_node(char *ensembleID, char *gene_name)
{
    if (ensembleID == NULL || gene_name == NULL)
    {
        return NULL;
    }

    geneID_name_node *new_node = (geneID_name_node *)malloc(sizeof(geneID_name_node));
    new_node->ensembleID = strdup(ensembleID);
    new_node->gene_name = strdup(gene_name);
    new_node->left = NULL;
    new_node->right = NULL;
    return new_node;
}

geneID_name_node *insert_geneID_name_node(
    geneID_name_node *root,
    char *ensembleID,
    char *gene_name)
{
    if (root == NULL)
    {
        return new_geneID_name_node(ensembleID, gene_name);
    }

    if (strcmp(ensembleID, root->ensembleID) < 0)
    {
        root->left = insert_geneID_name_node(root->left, ensembleID, gene_name);
    }

    if (strcmp(ensembleID, root->ensembleID) > 0)
    {
        root->right = insert_geneID_name_node(root->right, ensembleID, gene_name);
    }

    return root;
}

size_t count_geneID_name_node(geneID_name_node *root)
{
    if (root == NULL)
    {
        return 0;
    }

    return 1 + count_geneID_name_node(root->left) + count_geneID_name_node(root->right);
}

// get inorder traversal index of a geneID

void free_geneID_name_node(geneID_name_node *root)
{
    if (root == NULL)
    {
        return;
    }

    free_geneID_name_node(root->left);
    free_geneID_name_node(root->right);
    free(root->ensembleID);
    free(root->gene_name);
    free(root);
}

// assume node is sorted and ensembleID is in the node
void inorder_index_geneID_name_node(geneID_name_node *root, char *ensembleID, int *index)
{
    if (root == NULL)
    {
        return;
    }

    inorder_index_geneID_name_node(root->left, ensembleID, index);
    if (strcmp(root->ensembleID, ensembleID) <= 0)
    {
        (*index)++;
    }
    inorder_index_geneID_name_node(root->right, ensembleID, index);
}

int search_geneID_for_rownumber(geneID_name_node *root, char *ensembleID)
{
    int index = 0;
    inorder_index_geneID_name_node(root, ensembleID, &index);
    return index;
}

void write_gene_UMI_node_output(
    size_t cell_counts,
    gene_UMI_node *gene_UMI_root,
    geneID_name_node *geneID_name_root,
    gzFile fp)
{
    if (gene_UMI_root == NULL)
    {
        return;
    }

    write_gene_UMI_node_output(cell_counts, gene_UMI_root->left, geneID_name_root, fp);

    int rownumber = search_geneID_for_rownumber(geneID_name_root, gene_UMI_root->geneID);
    size_t unique_UMI_counts = count_UMI_node(gene_UMI_root->UMI);

    gzprintf(fp, "%d %zu %zu\n", rownumber, cell_counts, unique_UMI_counts);

    write_gene_UMI_node_output(cell_counts, gene_UMI_root->right, geneID_name_root, fp);
}

void write_cell_UMI_node_output(
    size_t *cell_counts,
    cell_gene_node *cell_gene_root,
    geneID_name_node *geneID_name_root,
    gzFile fp_barcodes,
    gzFile fp_matrix)
{

    if (cell_gene_root == NULL)
    {
        return;
    }

    write_cell_UMI_node_output(
        cell_counts,
        cell_gene_root->left,
        geneID_name_root,
        fp_barcodes,
        fp_matrix);

    // write barcodes
    gzprintf(fp_barcodes, "%s\n", cell_gene_root->CB);

    write_gene_UMI_node_output(
        *cell_counts,
        cell_gene_root->gene_UMI,
        geneID_name_root,
        fp_matrix);

    (*cell_counts)++;

    write_cell_UMI_node_output(
        cell_counts,
        cell_gene_root->right,
        geneID_name_root,
        fp_barcodes,
        fp_matrix);
}

void write_geneID_name_node_output(
    geneID_name_node *geneID_name_root,
    gzFile fp)
{
    if (geneID_name_root == NULL)
    {
        return;
    }

    write_geneID_name_node_output(geneID_name_root->left, fp);

    gzprintf(fp, "%s\t%s\tGene Expression\n", geneID_name_root->ensembleID, geneID_name_root->gene_name);

    write_geneID_name_node_output(geneID_name_root->right, fp);
}
