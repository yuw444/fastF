#include "count.h"

long total_n_features;
long total_cells;
long total_genes;

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

UMI_node *new_UMI_node(char *UMI)
{
    if (UMI == NULL)
    {
        return NULL;
    }
    UMI_node *root = (UMI_node *)malloc(sizeof(UMI_node));
    root->UMI = strdup(UMI);
    total_n_features++;
    root->left = NULL;
    root->right = NULL;
    return root;
}

UMI_node *insert_UMI_node(UMI_node *root, char *UMI)
{
    if (root == NULL)
    {
        root = new_UMI_node(UMI);
        return root;
    }
    else
    {
        int cmp = strcmp(UMI, root->UMI);
        if (cmp < 0)
        {
            root->left = insert_UMI_node(root->left, UMI);
        }
        if (cmp > 0)
        {
            root->right = insert_UMI_node(root->right, UMI);
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

    if (root->UMI != NULL)
    {
        free(root->UMI);
    }

    free_UMI_node(root->left);
    free_UMI_node(root->right);
    free(root);
}

gene_UMI_node *new_gene_UMI_node(size_t nthGene, char *UMI)
{
    if (UMI == NULL)
    {
        return NULL;
    }
    gene_UMI_node *root = (gene_UMI_node *)malloc(sizeof(gene_UMI_node));
    root->nthGene = nthGene;
    root->len_UMIs = 0;
    root->UMI = new_UMI_node(UMI);
    root->left = NULL;
    root->right = NULL;
    return root;
}

gene_UMI_node *insert_gene_UMI_node(
    gene_UMI_node *root,
    size_t nthGene,
    char *UMI)
{
    if (root == NULL)
    {
        root = new_gene_UMI_node(nthGene, UMI);
        return root;
    }
    else
    {
        int cmp = nthGene - root->nthGene;
        if (cmp < 0)
        {
            root->left = insert_gene_UMI_node(root->left, nthGene, UMI);
        }
        else if (cmp > 0)
        {
            root->right = insert_gene_UMI_node(root->right, nthGene, UMI);
        }
        else
        {
            root->UMI = insert_UMI_node(root->UMI, UMI);
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

    if (root->UMI != NULL)
    {
        free_UMI_node(root->UMI);
    }

    free_gene_UMI_node(root->left);
    free_gene_UMI_node(root->right);

    free(root);
}

cell_gene_node *new_cell_gene_node(size_t nthCell, size_t nthGene, char *UMI)
{
    if (UMI == NULL)
    {
        return NULL;
    }
    cell_gene_node *root = (cell_gene_node *)malloc(sizeof(cell_gene_node));
    root->nthCell = nthCell;
    root->gene_UMI = new_gene_UMI_node(nthGene, UMI);
    root->left = NULL;
    root->right = NULL;
    return root;
}

cell_gene_node *insert_cell_gene_node(cell_gene_node *root, size_t nthCell, size_t nthGene, char *UMI)
{
    if (root == NULL)
    {
        root = new_cell_gene_node(nthCell, nthGene, UMI);
        return root;
    }
    else
    {
        size_t cmp = nthCell - root->nthCell;
        if (cmp < 0)
        {
            root->left = insert_cell_gene_node(root->left, nthCell, nthGene, UMI);
        }
        else if (cmp > 0)
        {
            root->right = insert_cell_gene_node(root->right, nthCell, nthGene, UMI);
        }
        else
        {
            root->gene_UMI = insert_gene_UMI_node(root->gene_UMI, nthGene, UMI);
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

    if (root->gene_UMI != NULL)
    {
        free_gene_UMI_node(root->gene_UMI);
    }

    free_cell_gene_node(root->left);
    free_cell_gene_node(root->right);

    free(root);
}

list_saver_node *new_list_saver_node(char *element, size_t nthElement)
{
    if (element == NULL)
    {
        return NULL;
    }

    list_saver_node *root = (list_saver_node *)malloc(sizeof(list_saver_node));
    root->element = strdup(element);
    root->nthElement = nthElement;
    root->left = NULL;
    root->right = NULL;

    return root;
}

list_saver_node *insert_list_saver_node(list_saver_node *root, char *element, size_t nthElement)
{
    if (root == NULL)
    {
        root = new_list_saver_node(element, nthElement);
        return root;
    }
    else
    {
        int cmp = strcmp(element, root->element);
        if (cmp < 0)
        {
            root->left = insert_list_saver_node(root->left, element, nthElement);
        }
        else if (cmp > 0)
        {
            root->right = insert_list_saver_node(root->right, element, nthElement);
        }
        else
        {
            return root;
        }
    }
    return root;
}

size_t search_list_saver_node(list_saver_node *root, char *element, size_t len_elements)
{
    if (root == NULL)
    {
        return 0;
    }
    else
    {
        int cmp = strncmp(element, root->element, len_elements);
        if (cmp < 0)
        {
            return search_list_saver_node(root->left, element, len_elements);
        }
        else if (cmp > 0)
        {
            return search_list_saver_node(root->right, element, len_elements);
        }
        else
        {
            return root->nthElement;
        }
    }
}

void free_list_saver_node(list_saver_node *root)
{
    if (root == NULL)
    {
        return;
    }

    free(root->element);

    free_list_saver_node(root->left);
    free_list_saver_node(root->right);

    free(root);
}

cell_gene_node *sample_bam_UMI(
    char *bam_file,
    char *feature_file,
    char *barcode_file,
    double rate_reads,
    unsigned int seed)
{
    srand(seed);

    // initialize barcode tree
    node *barcode_tree = NULL;

    FILE *fp_barcode = fopen(barcode_file, "r");
    if (fp_barcode == NULL)
    {
        printf("ERROR: Cannot open Barcode file %s\n", barcode_file);
        exit(1);
    }
    char CB[MAX_LINE_LENGTH];
    while (fgets(CB, MAX_LINE_LENGTH, fp_barcode) != NULL)
    {
        total_cells++;
        // remove newline character
        CB[16] = '\0';
        barcode_tree = insert_tree(barcode_tree, CB);
    fclose(fp_barcode);

    // write the order of barcodes to a file
    gzFile fp_barcode_order = gzopen("barcode.tsv.gz", "w"); 
    if (fp_barcode_order == NULL)
    {
        printf("ERROR: Cannot open Barcode file %s\n", "barcode.tsv.gz");
        exit(1);
    }
    print_tree_gz(barcode_tree, fp_barcode_order);
    gzclose(fp_barcode_order);

    // initialize gene ID tree
    list_saver_node *geneID_tree = NULL;
    FILE *fp_feature = fopen(feature_file, "r");
    if (fp_feature == NULL)
    {
        printf("ERROR: Cannot open Feature file %s\n", feature_file);
        exit(1);
    }
    char line[MAX_LINE_LENGTH];
    while (fgets(line, MAX_LINE_LENGTH, fp_feature) != NULL)
    {
        total_genes++;
        char *geneID = strtok(line, "\t");
        // printf("geneID: %s\n", geneID);
        geneID_tree = insert_list_saver_node(geneID_tree, geneID, total_genes);
    }
    fclose(fp_feature);


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

        double random_number = (double)rand() / (double)RAND_MAX;

        // as some read may not be aligned to any gene, we need to check if TX tag exists
        if (cb != NULL)
        {
            char *cell_barcode = bam_aux2Z(cb);
            // printf("%s\n", cell_barcode);

            if (in(barcode_tree, cell_barcode, 16))
            {
                if (random_number < rate_reads)
                {
                    valid_read_counts++;

                    uint8_t *xf = bam_aux_get(bam_record, "xf");
                    int UMI_quality = bam_aux2i(xf);

                    if (UMI_quality == 25)
                    {
                        uint8_t *gx = bam_aux_get(bam_record, "GX");
                        uint8_t *ub = bam_aux_get(bam_record, "UB");
                        char *gene_ID = bam_aux2Z(gx);
                        char *UMI = bam_aux2Z(ub);

                        int nthCell = 0;
                        get_node_inorder_index(barcode_tree, cell_barcode, 16, &nthCell);

                        size_t nthGene = search_list_saver_node(geneID_tree, gene_ID, 18);

                        if (nthCell != 0 && nthGene != 0)
                        {
                            cell_gene_root = insert_cell_gene_node(
                                cell_gene_root,
                                nthCell,
                                nthGene,
                                UMI);
                        }
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

    if (barcode_tree != NULL)
    {
        free_list_saver_node(barcode_tree);
    }

    if (geneID_tree != NULL)
    {
        free_list_saver_node(geneID_tree);
    }

    return cell_gene_root;
}

void write_gene_UMI_node_output(
    size_t nthCell,
    gene_UMI_node *gene_UMI_root,
    gzFile fp)
{
    if (gene_UMI_root == NULL)
    {
        return;
    }

    write_gene_UMI_node_output(nthCell, gene_UMI_root->left, fp);

    gzprintf(fp, "%zu %zu %zu\n", gene_UMI_root->nthGene, nthCell, gene_UMI_root->len_UMIs);

    write_gene_UMI_node_output(nthCell, gene_UMI_root->right, fp);
}

void count_



void write_cell_gene_node_output(
    cell_gene_node *cell_gene_root,
    gzFile fp_matrix)
{

    if (cell_gene_root == NULL)
    {
        return;
    }

    write_cell_gene_node_output(
        cell_gene_root->left,
        fp_matrix);

    write_gene_UMI_node_output(
        cell_gene_root->nthCell,
        cell_gene_root->gene_UMI,
        fp_matrix);

    write_cell_gene_node_output(
        cell_gene_root->right,
        fp_matrix);
}

bool in_list(list_saver_node *root, char *element, int len_elements)
{
    if (root == NULL)
    {
        return 0;
    }

    if (strncmp(element, root->element, len_elements) == 0)
    {
        return 1;
    }
    else if (strncmp(element, root->element, len_elements) < 0)
    {
        return in_list(root->left, element, len_elements);
    }
    else
    {
        return in_list(root->right, element, len_elements);
    }
}

void print_list_saver_node(list_saver_node *root)
{
    if (root == NULL)
    {
        return;
    }

    print_list_saver_node(root->left);
    printf("%s\n", root->element);
    print_list_saver_node(root->right);
}