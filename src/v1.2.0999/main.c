#include <dirent.h>
#include <sys/stat.h>
#include "argparse.h"
#include "filter.h"
#include "count.h"
#include "extract.h"
#include "bam2db_ds.h"

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

int cmd_freq(int argc, const char **argv)
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
                      "\nFilter fastq file using cell barcode whitelist and read depth.",
                      "");

    argc = argparse_parse(&argparse, argc, argv);

    printf("whitelist: %s\n", whitelist_arg);

    if (path_R1_arg == NULL)
    {
        printf("\x1b[31mError:\x1b[0m path to R1 fastq files can not been NULL while filtering .\n");
        exit(1);
    }

    if (whitelist_arg == NULL && !all_cell)
    {
        printf("\x1b[31mError:\x1b[0m whitelist and --all cell option can not been both NULL at the same time.\n");
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

int cmd_crb(int argc, const char **argv)
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
        "\nExtract CR and CB tags from bam file and summarize them to a tsv file.",
        ""
    );

    argc = argparse_parse(&argparse, argc, argv);

    if (path_bam_arg == NULL)
    {
        printf("\x1b[31mError:\x1b[0m path to bam file can not been NULL while extracting .\n");
        exit(1);
    }

    char *path_out = (char *)malloc(1024 * sizeof(char));
    sprintf(path_out, "%s/CR_CB.tsv.gz", path_out_arg);

    // open file stream
    gzFile file_out = gzopen(path_out, "w");
    if (file_out == NULL)
    {
        printf("\x1b[31mError:\x1b[0m can not open file %s\n", path_out);
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

int cmd_bam2db(int argc, const char **argv)
{
    char *path_bam_arg = NULL;
    char *path_feature_arg = NULL;
    char *path_barcode_arg = NULL;
    float rate_cell;
    float rate_depth;
    char *name_database_arg = NULL;
    char *path_out_arg = ".";
    unsigned int seed_arg = 926;

    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_STRING('b', "bam", &path_bam_arg, "path to bam file", NULL, 0, 0),
        OPT_STRING('f', "feature", &path_feature_arg, "path to feature list file", NULL, 0, 0),
        OPT_STRING('a', "barcode", &path_barcode_arg, "path to barcode list file", NULL, 0, 0),
        OPT_STRING('d', "dbname", &name_database_arg, "name of database", NULL, 0, 0),
        OPT_FLOAT('c', "cell", &rate_cell, "rate of cell barcode", NULL, 0, 0),
        OPT_FLOAT('r', "depth", &rate_depth, "rate of depth", NULL, 0, 0),
        OPT_STRING('o', "out", &path_out_arg, "path to output directory", NULL, 0, 0),
        OPT_INTEGER('s', "seed", &seed_arg, "seed for random number generator", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usages, 0);
    argparse_describe(
        &argparse,
        "\nConvert bam files to SQL database for further downsample query.",
        ""
    );

    argc = argparse_parse(&argparse, argc, argv);

    // check the existence of input files
    if (access(path_bam_arg, F_OK) == -1)
    {
        printf("\x1b[31mError:\x1b[0m bam file: %s does not exist.\n", path_bam_arg);
        exit(1);
    }

    if (access(path_feature_arg, F_OK) == -1)
    {
        printf("\x1b[31mError:\x1b[0m feature file: %s does not exist.\n", path_feature_arg);
        exit(1);
    }

    if (access(path_barcode_arg, F_OK) == -1)
    {
        printf("\x1b[31mError:\x1b[0m barcode file: %s does not exist.\n", path_barcode_arg);
        exit(1);
    }

    if (access(name_database_arg, F_OK) != -1)
    {
        printf("\x1b[31mError:\x1b[0m database: %s already exists, change the name of database in -d argument!\n", name_database_arg);
        exit(1);
    }

    bam2db(
        path_bam_arg, 
        name_database_arg, 
        path_barcode_arg, 
        path_feature_arg, 
        rate_cell,
        rate_depth,
        seed_arg);
    
    // open the database
    sqlite3 *db;
    int rc = sqlite3_open(name_database_arg, &db);

    if (rc)
    {
        printf("\x1b[31mError:\x1b[0m can not open database: %s\n", name_database_arg);
        return 1;
    }

    // file path of output
    char *path_barcode = (char *)malloc(1024 * sizeof(char));
    sprintf(path_barcode, "%s/barcodes.tsv.gz", path_out_arg);
    gzFile file_barcode = gzopen(path_barcode, "wb");
    if (file_barcode == NULL)
    {
        printf("\x1b[31mError:\x1b[0m can not open file %s\n", path_barcode);
        exit(1);
    }
    char *path_feature = (char *)malloc(1024 * sizeof(char));
    sprintf(path_feature, "%s/features.tsv.gz", path_out_arg);
    gzFile file_feature = gzopen(path_feature, "wb");
    if (file_feature == NULL)
    {
        printf("\x1b[31mError:\x1b[0m can not open file %s\n", path_feature);
        exit(1);
    }
    char *path_matrix = (char *)malloc(1024 * sizeof(char));
    sprintf(path_matrix, "%s/matrix.mtx.gz", path_out_arg);
    gzFile file_matrix = gzopen(path_matrix, "wb");
    if (file_matrix == NULL)
    {
        printf("\x1b[31mError:\x1b[0m can not open file %s\n", path_matrix);
        exit(1);
    }

    // create summary table of umi
    char *error_msg = NULL;
    char *sql =   "CREATE TABLE mtx AS "
                  "SELECT feature_index, cell_index, COUNT(DISTINCT encoded_umi) AS expression_level "
                  " FROM umi "
                  " GROUP BY feature_index, cell_index;"; 
    if (sqlite3_exec(db, sql, NULL, 0, &error_msg) != SQLITE_OK)
    {
        printf("\x1b[31mError:\x1b[0m SQL error: %s\n", error_msg);
        sqlite3_free(error_msg);
        return 1;
    }

    // get the number of rows of every table
    size_t nrow_mtx = nrow_sql_table(db, "mtx");
    size_t nrow_barcode = nrow_sql_table(db, "cell");
    size_t nrow_feature = nrow_sql_table(db, "feature");
    
    // write the nrows to the header of matrix.mtx.gz
    char *header = (char *)malloc(1024 * sizeof(char));
    gzprintf(file_matrix, "%%MatrixMarket matrix coordinate integer general\n");
    gzprintf(file_matrix, "%%metadata_json: {\"software_version\": \"bamDesign-1.0.0\", \"format_version\": 1}\n");
    gzprintf(file_matrix, "%zu %zu %zu\n", nrow_feature, nrow_barcode, nrow_mtx);

    // write table mtx to matrix.mtx.gz
    table2gz(db, "mtx", file_matrix, 0, " ");
    printf("matrix.mtx.gz is generated.\n");

    // write table cell to barcodes.tsv.gz
    table2gz(db, "cell", file_barcode, 0, " ");
    printf("barcodes.tsv.gz is generated.\n");

    // write table feature to features.tsv.gz
    table2gz(db, "feature", file_feature, 0, "\t");
    printf("features.tsv.gz is generated.\n");

    free(path_barcode);
    free(path_feature);
    free(path_matrix);
    free(header);
    sqlite3_close(db);

    return 0;


}


static struct cmd_struct commands[] = {
    {"filter", cmd_filter},
    {"freq", cmd_freq},
    {"crb", cmd_crb},
    {"bam2db", cmd_bam2db},
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
