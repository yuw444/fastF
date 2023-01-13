#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "argparse.h"

static const char *const usages[] = {
    "fastF <OPTION>",
    NULL,
    NULL,
};

#define PERM_READ (1 << 0)
#define PERM_WRITE (1 << 1)
#define PERM_EXEC (1 << 2)

int main(int argc, const char **argv)
{
    const char *whitelist_arg = NULL;
    float rate_arg = 0.f;
    const char *path_i_arg = ".";
    const char *path_o_arg = ".";
    const char *sample_arg = NULL;
    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_GROUP("Basic options"),
        OPT_STRING('w', "whitelist", &whitelist_arg, "whitelist of cell barcodes", NULL, 0, 0),
        OPT_FLOAT('r', "rate", &rate_arg, "rate of reads to keep after matching cell_barcode", NULL, 0, 0),
        OPT_STRING('i', "in", &path_i_arg, "dir to sample fastq files", NULL, 0, 0),
        OPT_STRING('o', "out", &path_o_arg, "dir to output fastq files", NULL, 0, 0),
        OPT_STRING('s', "sample", &sample_arg, "sample name of fastq files", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usages, 0);
    argparse_describe(&argparse,
                      "\nFilter fastq file using cell barcode whitelist.",
                      "");

    argc = argparse_parse(&argparse, argc, argv);

    char **file_in_arg = malloc(3 * sizeof(char *));

    int i = 0;

    if (whitelist_arg != NULL && rate_arg != 0 && sample_arg != NULL)
    {
        DIR *dir_i, *dir_o;
        struct dirent *ptr;

        dir_i = opendir(path_i_arg);
        dir_o = opendir(path_o_arg);

        if (dir_i == NULL || dir_o == NULL)
        {
            printf("Error: cannot open dir %s or %s \n", path_i_arg, path_o_arg);
            exit(1);
        }
        else
        {
            while ((ptr = readdir(dir_i)) != NULL)
            {
                if (ptr->d_type == DT_REG)
                {
                    if (strncmp(ptr->d_name, sample_arg, strlen(sample_arg)) == 0 &&
                        strncmp(ptr->d_name + strlen(ptr->d_name) - 9, ".fastq.gz", 9) == 0)
                    {
                        file_in_arg[i] = malloc(strlen(ptr->d_name) + 1);
                        strcpy(file_in_arg[i], ptr->d_name);
                        i++;
                        // printf("%s \n", file_in_arg[i - 1]);
                    }

                    // printf("i is %d \n", i);

                    if (i > 3)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            free(file_in_arg[j]);
                        }
                        free(file_in_arg);
                        printf("Error: more than 3 fastq files with name start with %s in %s \n", sample_arg, path_i_arg);
                        closedir(dir_i);
                        closedir(dir_o);
                        exit(1);
                    }
                }
            }

            closedir(dir_i);
            closedir(dir_o);
        }
    }
    else
    {
        printf("Error: missing arguments for whitelist, rate or sample\n");
        exit(1);
    }

    char *I1, *R1, *R2;

    for(int j = 0; j < 3; j++)
    {
        if(strstr(file_in_arg[j], "R1") != NULL)
        {
            R1  = file_in_arg[j];
        }
        else if(strstr(file_in_arg[j], "R2") != NULL)
        {
            R2 = file_in_arg[j];
        }
        else if(strstr(file_in_arg[j], "I1") != NULL)
        {
            I1 = file_in_arg[j];
        }
        else
        {
            printf("Error: fastq file name is not in the format of sample_*_I1_*_.fastq.gz, sample_*_R1_*_.fastq.gz or sample_*_R2_*_.fastq.gz \n");
            exit(1);
        }
    }


    printf("whitelist is %s \n", whitelist_arg);
    printf("rate is %f \n", rate_arg);
    printf("path_i is %s \n", path_i_arg);
    printf("path_o is %s \n", path_o_arg);
    printf("sample is %s \n", sample_arg);
    printf("I1 is %s \n", I1);
    printf("R1 is %s \n", R1);
    printf("R2 is %s \n", R2);
    

    for (int j = 0; j < 3; j++)
    {
        free(file_in_arg[j]);
    }
    free(file_in_arg);

    // for (int j = 0; j < 3; j++)
    // {
    //     printf("%s\n", file_in_arg[j]);
    // }

    return 0;
}