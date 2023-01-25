#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include "argparse.h"
#include "fastq_filter.h"

static const char *const usages[] = {
    "fastF <OPTION>",
    NULL,
    NULL,
};

int main(int argc, const char **argv)
{
    double startTime = clock();

    char *whitelist_arg = NULL;
    int seed_arg = 926;
    float rate_arg = 0.f;
    const char *path_i_arg = ".";
    const char *path_o_arg = ".";
    const char *sample_arg = NULL;
    struct argparse_option options[] = {
        OPT_HELP(),
        OPT_GROUP("Basic options"),
        OPT_STRING('w', "whitelist", &whitelist_arg, "whitelist of cell barcodes", NULL, 0, 0),
        OPT_INTEGER('s', "seed", &seed_arg, "seed for random number generator", NULL, 0, 0),
        OPT_FLOAT('r', "rate", &rate_arg, "rate of reads to keep after matching cell_barcode", NULL, 0, 0),
        OPT_STRING('i', "in", &path_i_arg, "dir to sample fastq files", NULL, 0, 0),
        OPT_STRING('o', "out", &path_o_arg, "dir to output fastq files", NULL, 0, 0),
        OPT_STRING('n', "sample", &sample_arg, "sample name of fastq files", NULL, 0, 0),
        OPT_END(),
    };

    struct argparse argparse;
    argparse_init(&argparse, options, usages, 0);
    argparse_describe(&argparse,
                      "\nFilter fastq file using cell barcode whitelist.",
                      "");

    argc = argparse_parse(&argparse, argc, argv);

    char **file_in_names = malloc(3 * sizeof(char *));

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
            char *output_fold = malloc(strlen(path_o_arg) + strlen(sample_arg) + 2);
            sprintf(output_fold, "%s/%s", path_o_arg, sample_arg);

            errno = 0;
            int ret = mkdir(output_fold, 0777);
            free(output_fold);
            if (ret == -1)
            {
                switch (errno)
                {
                case EACCES:
                    printf("Error: The parent directory %s does not allow write\n", path_o_arg);
                    exit(EXIT_FAILURE);
                case EEXIST:
                    printf("Error: Folder %s already exists in %s.\n", sample_arg, path_o_arg);
                    exit(EXIT_FAILURE);
                case ENAMETOOLONG:
                    printf("Error: Pathname %s is too long\n", sample_arg);
                    exit(EXIT_FAILURE);
                default:
                    perror("mkdir\n");
                    exit(EXIT_FAILURE);
                }
            }

            while ((ptr = readdir(dir_i)) != NULL)
            {
                if (ptr->d_type == DT_REG)
                {
                    if (strncmp(ptr->d_name, sample_arg, strlen(sample_arg)) == 0 &&
                        strncmp(ptr->d_name + strlen(ptr->d_name) - 9, ".fastq.gz", 9) == 0)
                    {
                        file_in_names[i] = malloc(strlen(ptr->d_name) + 1);
                        strcpy(file_in_names[i], ptr->d_name);
                        i++;
                        // printf("%s \n", file_in_names[i - 1]);
                    }

                    // printf("i is %d \n", i);
                }
            }

            if (i > 3)
            {
                for (int j = 0; j < i; j++)
                {
                    free(file_in_names[j]);
                }
                free(file_in_names);

                printf("Error: more than 3 fastq files with name start with %s in %s \n", sample_arg, path_i_arg);
                closedir(dir_i);
                closedir(dir_o);
                exit(1);
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

    for (int j = 0; j < 3; j++)
    {
        if (strstr(file_in_names[j], "R1") != NULL)
        {
            R1 = file_in_names[j];
        }
        else if (strstr(file_in_names[j], "R2") != NULL)
        {
            R2 = file_in_names[j];
        }
        else if (strstr(file_in_names[j], "I1") != NULL)
        {
            I1 = file_in_names[j];
        }
        else
        {
            printf("Error: fastq file name is not in the format of sample_*_I1_*_.fastq.gz, sample_*_R1_*_.fastq.gz or sample_*_R2_*_.fastq.gz \n");
            exit(1);
        }
    }

    char **file_in_path = malloc(3 * sizeof(char *));
    char **file_out_path = malloc(3 * sizeof(char *));

    for (int j = 0; j < 3; j++)
    {
        file_in_path[j] = malloc(strlen(path_i_arg) + strlen(file_in_names[j]) + 100);
        file_out_path[j] = malloc(strlen(path_o_arg) + strlen(file_in_names[j]) + 100);
    }

    sprintf(file_in_path[0], "%s/%s", path_i_arg, I1);
    sprintf(file_in_path[1], "%s/%s", path_i_arg, R1);
    sprintf(file_in_path[2], "%s/%s", path_i_arg, R2);

    sprintf(file_out_path[0], "%s/%s/%s", path_o_arg, sample_arg, I1);
    sprintf(file_out_path[1], "%s/%s/%s", path_o_arg, sample_arg, R1);
    sprintf(file_out_path[2], "%s/%s/%s", path_o_arg, sample_arg, R2);

    printf("whitelist is %s \n", whitelist_arg);
    printf("rate is %f \n", rate_arg);
    printf("path_i is %s \n", path_i_arg);
    printf("path_o is %s \n", path_o_arg);
    printf("sample is %s \n", sample_arg);
    printf("I1 is %s \n", file_in_path[0]);
    printf("R1 is %s \n", file_in_path[1]);
    printf("R2 is %s \n", file_in_path[2]);
    printf("I1_out is %s \n", file_out_path[0]);
    printf("R1_out is %s \n", file_out_path[1]);
    printf("R2_out is %s \n", file_out_path[2]);

    /********************************/
    /*      process fastq files     */
    /********************************/

    for (int j = 0; j < 3; j++)
    {
        file_in[j] = gzopen(file_in_path[j], "rb");
        gzbuffer(file_in[j], 128 * 1024);
        file_out[j] = gzopen(file_out_path[j], "wb");
        gzbuffer(file_out[j], 128 * 1024);
    }

    reader_buffer = init_queue();
    writer_buffer = init_queue();

    rate_threshold = rate_arg;

    printf("Reading whitelist...\n");
    int nrow = get_row(whitelist_arg);
    printf("nrow = %d\n", nrow);
    char **whitelist = read_txt(whitelist_arg, nrow);

    qsort(whitelist, nrow, sizeof(char *), compare);
    tree_whitelist = construct_tree(whitelist, 0, nrow - 1);

    // print_tree(tree_whitelist);
    int t1 = in(tree_whitelist, "TTTGGTTGTGACCAAG");
    printf("t1 = %d\n", t1);

    int *reader_id = malloc(NUM_READERS * sizeof(int));
    int *proc_id = malloc(NUM_PROCESSORS * sizeof(int));
    int *writer_id = malloc(NUM_WRITERS * sizeof(int));

    srand(seed_arg);

    pthread_t reader_thread[NUM_READERS];
    pthread_t processor_thread[NUM_PROCESSORS];
    pthread_t writer_thread[NUM_WRITERS];

    for (int i = 0; i < NUM_READERS; i++)
    {
        reader_id[i] = i;
        pthread_create(&reader_thread[i], NULL, reader, &reader_id[i]);
    }

    sleep(1);

    for (int i = 0; i < NUM_PROCESSORS; i++)
    {
        proc_id[i] = i;
        pthread_create(&processor_thread[i], NULL, processor, &proc_id[i]);
    }

    for (int i = 0; i < NUM_WRITERS; i++)
    {
        writer_id[i] = i;
        pthread_create(&writer_thread[i], NULL, writer, &writer_id[i]);
    }

    // wait for the reader and processor threads to finish
    for (int i = 0; i < NUM_READERS; i++)
    {
        pthread_join(reader_thread[i], NULL);
    }

    for (int i = 0; i < NUM_PROCESSORS; i++)
    {
        pthread_join(processor_thread[i], NULL);
    }

    for (int i = 0; i < NUM_WRITERS; i++)
    {
        pthread_join(writer_thread[i], NULL);
    }

    free(proc_id);
    free(reader_id);
    free(writer_id);
    free(reader_buffer);
    free(writer_buffer);

    free_tree_node(tree_whitelist);
    free(whitelist);

    for (int i = 0; i < 3; i++)
    {
        gzclose(file_in[i]);
        gzclose(file_out[i]);
    }

    for (int j = 0; j < 3; j++)
    {
        free(file_in_names[j]);
        free(file_in_path[j]);
        free(file_out_path[j]);
    }

    free(file_in_names);
    free(file_in_path);
    free(file_out_path);

    double stopTime = clock();

    double secsElapsed = (stopTime - startTime) / CLOCKS_PER_SEC;

    printf("Elapsed: %f seconds\n", secsElapsed);
    return 0;
}

// ./test_main -w ../data/whitelist.txt -r 1.0 -i ../data -o ../data/ -n test

// ./test_main -w /home/rstudio/Frag/data/whitelist.txt -r 1.0 -i /home/rstudio/Frag/data -o /home/rstudio/Frag/data/ -n test

// Elapsed: 397.7s multi-threading with 4 threads
// Elapsed: 307.436405 seconds used `gzbuffer` to speed up with 1 consumer and 128k buffer
// Elapsed: 336.436405 seconds used `gzbuffer` to speed up with 2 consumer and 128k buffer
// Elapsed: 378.906889 seconds used `gzbuffer` to speed up with 4 consumer and 128k buffer
// Elapsed: 275.436405 seconds used `gzbuffer` to speed up with 1 reader, 1 processor, 1 writer and 128k buffer
// Elapsed: 374.329000 seconds used `gzbuffer` to speed up with 1 reader, 3 processor, 1 writer and 128k buffer
// Elapsed: 374.329000 seconds used `gzbuffer` to speed up with 1 reader, 3 processor, 2 writer and 128k buffer