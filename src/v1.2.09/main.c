#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include "argparse.h"
#include "fastq_filter.h"
#include "count.h"

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

    if (path_R1_arg == NULL)
    {
        printf("Error: path to R1 fastq files can not been NULL while filtering .\n");
        exit(1);
    }

    if (whitelist_arg == NULL || !all_cell)
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

    file_in[0] = gzopen(path_I1_arg, "r");
    file_in[1] = gzopen(path_R1_arg, "r");
    file_in[2] = gzopen(path_R2_arg, "r");

    if (path_I1_arg == NULL)
    {
        file_in[0] = Z_NULL;
        file_out[0] = Z_NULL;
    }
    else
    {
       file_in[0] = gzopen(path_I1_arg, "r");
       file_out[0] = gzopen(strcat(path_o_arg, "/I1.fastq.gz"), "w"); 
    }

    file_in[1] = gzopen(path_R1_arg, "r");
    file_out[1] = gzopen(strcat(path_o_arg, "/R1.fastq.gz"), "w");

    if (path_R2_arg == NULL)
    {
        file_in[2] = Z_NULL;
        file_out[2] = Z_NULL;
    }
    else
    {
        file_in[2] = gzopen(path_R2_arg, "r");
        file_out[2] = gzopen(strcat(path_o_arg, "/R2.fastq.gz"), "w");
    }


    if (whitelist_arg != NULL)
    {
        printf("Reading whitelist...\n");
        int nrow = get_row(whitelist_arg);
        // printf("nrow = %d\n", nrow);
        char **whitelist = read_txt(whitelist_arg, nrow);

        tree_whitelist = construct_tree(whitelist, nrow);

        printf("Processing fastq files...\n");
        fastF(file_in, file_out, tree_whitelist, len_cellbarcode, seed_arg, rate_arg, all_cell);

        free_tree_node(tree_whitelist);
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
        gzclose(file_in[i]);
        gzclose(file_out[i]);
    }

    return 0;
}

static struct cmd_struct commands[] = {
    {"filter", cmd_filter},
    {"whitelist", cmd_whitelist},
};

int main(int argc, const char **argv)
{
    clock_t startTime = clock();

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

    clock_t endTime = clock();

    printf("Elapsed time: %f seconds\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);

    return 0;
}