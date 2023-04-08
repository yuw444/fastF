#include <dirent.h>
#include <sys/stat.h>
#include "argparse.h"
#include "fastq_filter.h"
#include "count.h"
#include "extract.h"

#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))

static const char *const usages[] = {
    "fastF subcommands <options>",
    NULL,
};

struct cmd_struct
{
    const char *cmd;
    int (*fn)(int, const char **);
};

int cmd_whitelist(int argc, const char **argv)
{

    char * path_R1_arg = NULL;
    char * path_out_arg = NULL;
    size_t len_cellbarcode = 16;
    size_t len_umi = 10;

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_GROUP("Basic options"),
        OPT_STRING('R', "R1", &path_R1_arg, "path to R1 fastq files", NULL, 0, 0),
        OPT_STRING('o', "out", &path_out_arg, "path to output whitelist", NULL, 0, 0),
        OPT_INTEGER('l', "len", &len_cellbarcode, "length of cell barcode", NULL, 0, 0),
        OPT_INTEGER('u', "umi", &len_umi, "length of UMI", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usages, 0);
    argparse_describe(&argparse,
                      "\nFind all the cell barcode whitelist and their frequencies.",
                      "");

    argc = argparse_parse(&argparse, argc, argv);
    
    if (path_R1_arg == NULL)
    {
        fprintf(stderr, "Please specify the path to R1 fastq files.\n");
        exit(1);
    }

    gzFile R1 = gzopen(path_R1_arg, "r");

    if (R1 == NULL)
    {
        fprintf(stderr, "Cannot open file %s \n", path_R1_arg);
        exit(1);
    }

    node *head = cell_counts(R1, len_cellbarcode, len_umi);

    gzclose(R1);

    // output the cell barcode whitelist
    char path_out[1024];
    sprintf(path_out, "%s/whitelist.txt", path_out_arg);
    FILE *fp = fopen(path_out, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Cannot open file %s \n", path_out);
        exit(1);
    }

    print_tree(head, fp);

    fclose(fp);

    // free the memory
    free_tree_node(head);

    return 0;

}

int cmd_filter(int argc, const char **argv)
{
    char *path_I1_arg = NULL;
    char *path_R1_arg = NULL;
    char *path_R2_arg = NULL;
    char *path_o_arg = ".";
    char *whitelist_arg = NULL;
    size_t len_cellbarcode = 16;
    int seed_arg = 926;
    float rate_arg = 0.f;
    int all_cell = 0;

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_GROUP("Basic options"),
        OPT_STRING('I', "I1", &path_I1_arg, "path to sample I1 fastq files", NULL, 0, 0),
        OPT_STRING('R', "R1", &path_R1_arg, "path to sample R1 fastq files", NULL, 0, 0),
        OPT_STRING('r', "R2", &path_R2_arg, "path to sample R2 fastq files", NULL, 0, 0),
        OPT_STRING('o', "out", &path_o_arg, "dir to output fastq files", NULL, 0, 0),
        OPT_STRING('w', "whitelist", &whitelist_arg, "whitelist of cell barcodes", NULL, 0, 0),
        OPT_INTEGER('l', "len", &len_cellbarcode, "length of cell barcode", NULL, 0, 0),
        OPT_INTEGER('s', "seed", &seed_arg, "seed for random number generator", NULL, 0, 0),
        OPT_FLOAT('t', "rate", &rate_arg, "rate of reads to keep after matching cell_barcode", NULL, 0, 0),
        OPT_BOOLEAN('a', "allcells", &all_cell, "keep all reads with cell barcode", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usages, 0);
    argparse_describe(&argparse,
                      "\nFilter fastq file using cell barcode whitelist.",
                      "");

    argc = argparse_parse(&argparse, argc, argv);

    printf("whitelist: %s\n", whitelist_arg);

    if (path_R1_arg == NULL)
    {
        printf("Error: path to R1 fastq files can not been NULL while filtering .\n");
        exit(1);
    }

    if (whitelist_arg == NULL && !all_cell)
    {
        printf("Error: whitelist and --all cell option can not been both NULL at the same time.\n");
        exit(1);
    }

    /********************************/
    /*      process fastq files     */
    /********************************/

    // shared data by the reader and processors
    //------------------------------------------------
    gzFile file_in[3];
    //------------------------------------------------

    // shared data by the processors and writer
    //------------------------------------------------
    gzFile file_out[3];
    node *tree_whitelist;
    //------------------------------------------------

    // file_in[0] = gzopen(path_I1_arg, "r");
    // file_in[1] = gzopen(path_R1_arg, "r");
    // file_in[2] = gzopen(path_R2_arg, "r");

    if (path_I1_arg == NULL)
    {
        file_in[0] = Z_NULL;
        file_out[0] = Z_NULL;
    }
    else
    {
       file_in[0] = gzopen(path_I1_arg, "r");
       char *path_I1_out = (char *)malloc(1024 * sizeof(char));
       sprintf(path_I1_out, "%s/I1.fastq.gz", path_o_arg);
       file_out[0] = gzopen(path_I1_out, "w"); 
       free(path_I1_out);
    }

    file_in[1] = gzopen(path_R1_arg, "r");
    char *path_R1_out = (char *)malloc(1024 * sizeof(char));
    sprintf(path_R1_out, "%s/R1.fastq.gz", path_o_arg);
    file_out[1] = gzopen(path_R1_out, "w");
    free(path_R1_out);

    if (path_R2_arg == NULL)
    {
        file_in[2] = Z_NULL;
        file_out[2] = Z_NULL;
        printf("TRUE\n");
    }
    else
    {
        file_in[2] = gzopen(path_R2_arg, "r");
        char *path_R2_out = (char *)malloc(1024 * sizeof(char));
        sprintf(path_R2_out, "%s/R2.fastq.gz", path_o_arg);
        file_out[2] = gzopen(path_R2_out, "w");
        free(path_R2_out);
    }

    if (whitelist_arg != NULL)
    {
        printf("Reading whitelist...\n");
        int nrow = get_row(whitelist_arg);
        printf("nrow = %d\n", nrow);
        char **whitelist = read_txt(whitelist_arg, nrow);

        tree_whitelist = construct_tree(whitelist, nrow);

        printf("Processing fastq files...\n");
        fastF(file_in, file_out, tree_whitelist, len_cellbarcode, seed_arg, rate_arg, all_cell);

        free_tree_node(tree_whitelist);
        for(int i = 0; i < nrow; i++)
        {
            free(whitelist[i]);
        }
        free(whitelist);
    } 
    else
    {   
        printf("Subsample fastq files directly without cell barcode whitelist...\n");
        printf("Processing fastq files...\n");
        fastF(file_in, file_out, NULL, len_cellbarcode, seed_arg, rate_arg, all_cell);
    }

    // print_tree(tree_whitelist);
    // int t1 = in(tree_whitelist, "TTTGGTTGTGACCAAG");
    // printf("t1 = %d\n", t1);

    for (int i = 0; i < 3; i++)
    {
        if(file_in[i] != Z_NULL)
        {
            gzclose(file_in[i]);
            gzclose(file_out[i]);
        }

    }

    return 0;
}

int cmd_extract(int argc, const char **argv)
{
    char *path_bam_arg = NULL;
    char *path_out_arg = ".";

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_STRING('b', "bam", &path_bam_arg, "path to bam file", NULL, 0, 0),
        OPT_STRING('o', "out", &path_out_arg, "path to output directory", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usages, 0);
    argparse_describe(
        &argparse,
        "\nExtract CR and CB from bam file.",
        ""
    );

    argc = argparse_parse(&argparse, argc, argv);

    if (path_bam_arg == NULL)
    {
        printf("Error: path to bam file can not been NULL while extracting .\n");
        exit(1);
    }

    char *path_out = (char *)malloc(1024 * sizeof(char));
    sprintf(path_out, "%s/CR_CB.tsv.gz", path_out_arg);

    // open file stream
    gzFile file_out = gzopen(path_out, "w");
    if (file_out == NULL)
    {
        printf("Error: can not open file %s\n", path_out);
        exit(1);
    }

    // read bam
    struct CB_node *CB_tree = read_bam(path_bam_arg);

    printf("Writing to file...\n");

    // write to file
    print_CB_node(CB_tree, file_out);

    // close file stream
    gzclose(file_out);

    // free memory
    free_CB_node(CB_tree);
    free(path_out);

    printf("Done.\n");
    return 0;
    
}

int cmd_sample(int argc, const char **argv)
{
    char *path_bam_arg = NULL;
    char *path_CB_arg = NULL;
    float rate_reads_arg = 1.0f;
    char *path_out_arg = ".";
    unsigned int all_cell = 0;
    int seed_arg = 926;

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_STRING('b', "bam", &path_bam_arg, "path to bam file", NULL, 0, 0),
        OPT_STRING('c', "CB", &path_CB_arg, "path to CB list file", NULL, 0, 0),
        OPT_FLOAT('r', "rate", &rate_reads_arg, "rate of reads to be sampled", NULL, 0, 0),
        OPT_STRING('o', "out", &path_out_arg, "path to output directory", NULL, 0, 0),
        OPT_BOOLEAN('a', "all", &all_cell, "sample all cells", NULL, 0, 0),
        OPT_INTEGER('s', "seed", &seed_arg, "seed for random number generator", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usages, 0);
    argparse_describe(
        &argparse,
        "\nExtract CR and CB from bam file.",
        ""
    );

    argc = argparse_parse(&argparse, argc, argv);

    if (path_bam_arg == NULL)
    {
        printf("Error: path to bam file can not been NULL while sampling .\n");
        exit(1);
    }

    if (path_CB_arg == NULL && all_cell == 0)
    {
        printf("Error: path to CB list file can not been NULL without all_cell option while sampling .\n");
        exit(1);
    }

    // open gzFile stream
    char *path_barcode = (char *)malloc(1024 * sizeof(char));
    char *path_feature = (char *)malloc(1024 * sizeof(char));
    char *path_matrix = (char *)malloc(1024 * sizeof(char));

    sprintf(path_barcode, "%s/barcodes.tsv.gz", path_out_arg);
    sprintf(path_feature, "%s/features.tsv.gz", path_out_arg);
    sprintf(path_matrix, "%s/matrix.mtx.gz", path_out_arg);

    
    gzFile file_barcode = gzopen(path_barcode, "w");
    gzFile file_feature = gzopen(path_feature, "w");
    gzFile file_matrix = gzopen(path_matrix, "w");

    if (file_barcode == Z_NULL || file_feature == Z_NULL || file_matrix == Z_NULL)
    {
        printf("Error: can not open output file in %s\n", path_out_arg);
        exit(1);
    }

    // read bam
    cell_gene_node_geneID_name_node *root = sample_bam_UMI(
        path_bam_arg,
        path_CB_arg, 
        all_cell,
        (double)rate_reads_arg,
        seed_arg);

    // count number of cells
    size_t n_cells = count_cell_gene_node(root->cell_gene_root);
    // count number of (unique) genes
    size_t n_genes = count_geneID_name_node(root->geneID_name_root);
    // count total number of non-zero genes across all cells
    size_t n_features = count_total_UMI_node(root->cell_gene_root);

    printf("Writing to file...\n");

    // write to matrix.mtx.gz
    gzprintf(file_matrix, "%%%%MatrixMarket matrix coordinate integer general\n");
    // gzprintf(file_matrix, "%%metadata_json: {\"software_version\": \"cellranger\", \"format_version\": 2}\n");
    gzprintf(file_matrix, "%%%%generated_by: fastF sample (v0.0.9)\n");
    gzprintf(file_matrix, "%zu %zu %zu\n", n_genes, n_cells, n_features);

    write_geneID_name_node_output(root->geneID_name_root, file_feature);
    size_t nth_cell = 1;

    write_cell_UMI_node_output(
        &nth_cell, 
        root->cell_gene_root, 
        root->geneID_name_root,
        file_barcode, 
        file_matrix
    );

    // close file stream
    gzclose(file_barcode);
    gzclose(file_feature);
    gzclose(file_matrix);

    // free memory
    free_cell_gene_node(root->cell_gene_root);
    free_geneID_name_node(root->geneID_name_root);
    free(root);
    free(path_barcode);
    free(path_feature);
    free(path_matrix);

    printf("Done.\n");
    return 0;
    
    // run command
    // ./fastF sample -b ../data/test.bam  -c whitelist.txt -a

}


static struct cmd_struct commands[] = {
    {"filter", cmd_filter},
    {"whitelist", cmd_whitelist},
    {"extract", cmd_extract},
    {"sample", cmd_sample},
};

int main(int argc, const char **argv)
{
    struct argparse argparse;
    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_END(),
    };

    argparse_init(&argparse, options, usages, ARGPARSE_STOP_AT_NON_OPTION);
    argc = argparse_parse(&argparse, argc, argv);
    if (argc < 1)
    {
        argparse_usage(&argparse);
        return -1;
    }

    /* Try to run command with args provided. */
    struct cmd_struct *cmd = NULL;
    for (int i = 0; i < ARRAY_SIZE(commands); i++)
    {
        if (!strcmp(commands[i].cmd, argv[0]))
        {
            cmd = &commands[i];
        }
    }

    if (cmd)
    {
        return cmd->fn(argc, argv);
    }

    return 0;
}
