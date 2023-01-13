// Author : Yu Wang
// Date : 01/06/2023
/*

This program is used to filter the file files based on the whitelist.
It is designed to be run with multiple threads.

*/

#include "fastq_filter.h"

#define NUM_CONSUMERS 4
#define NUM_PRODUCERS 1

// shared data by the reader and processors
//------------------------------------------------
queue *reader_buffer;
gzFile file_in[3];
//------------------------------------------------

// shared data by the processors and writer
//------------------------------------------------
pthread_mutex_t write_mutex = PTHREAD_MUTEX_INITIALIZER;
gzFile file_out[3];
node *tree_whitelist;
double rate_threshold;
//------------------------------------------------

// Function reader
void *producer(void *id)
{
    int producer_id = *((int *)id);
    printf("Producer thread %d is running!\n", producer_id);
    
    while (1)
    {
        while (!is_full(reader_buffer))
        {
            fastq *I1 = get_fastq(file_in[0]);
            if (flag == 1)
            {
                printf("Producer thread %d is exiting!\n", producer_id);
                pthread_exit(NULL);
            }
            fastq *R1 = get_fastq(file_in[1]);
            fastq *R2 = get_fastq(file_in[2]);

            // printf("Producer:\n %s%s%s", R1->id, R1->seq, R1->qual);

            comb_fastq *comb = malloc(sizeof(comb_fastq));
            comb->I1 = I1;
            comb->R1 = R1;
            comb->R2 = R2;

            enqueue(reader_buffer, comb);
        }
    }
}

void *consumer(void *id)
{
    int consumer_id = *((int *)id);
    printf("Consumer thread %d is running!\n", consumer_id);
    while (1)
    {
        while (!is_empty(reader_buffer))
        {
            comb_fastq *comb = dequeue(reader_buffer);

            fastq *I1 = comb->I1;
            fastq *R1 = comb->R1;
            fastq *R2 = comb->R2;

            // printf("Consumer:\n%s%s%s", comb->R1->id, comb->R1->seq, comb->R1->qual);

            char *cell_barcode = substring(R1->seq, 0, 16);

            if (in(tree_whitelist, cell_barcode))
            {
                pthread_mutex_lock(&write_mutex);
                char *comb_I1 = combine_string(I1);
                char *comb_R1 = combine_string(R1);
                char *comb_R2 = combine_string(R2);
                gzputs(file_out[0], comb_I1);
                gzputs(file_out[1], comb_R1);
                gzputs(file_out[2], comb_R2);
                free(comb_I1);
                free(comb_R1);
                free(comb_R2);
                pthread_mutex_unlock(&write_mutex);
            }
            free(cell_barcode);
            free_comb_fastq(comb);
        }
        if (flag == 1 && is_empty(reader_buffer))
        {
            printf("Consumer thread %d is exiting!\n", consumer_id);
            pthread_exit(NULL);
        }
    }

    // printf("processor is working!\n");
}

int main()
{
    file_in[0] = gzopen("../data/I1.fastq.gz", "rb");
    file_in[1] = gzopen("../data/R1.fastq.gz", "rb");
    file_in[2] = gzopen("../data/R2.fastq.gz", "rb");

    file_out[0] = gzopen("../data/filtered_I1.fastq.gz", "wb");
    file_out[1] = gzopen("../data/filtered_R1.fastq.gz", "wb");
    file_out[2] = gzopen("../data/filtered_R2.fastq.gz", "wb");

    reader_buffer = init_queue();

    rate_threshold = 1.0;

    char *filename_whitelist = "../data/whitelist.txt";

    printf("Reading whitelist...\n");
    int nrow = get_row(filename_whitelist);

    printf("nrow = %d\n", nrow);

    char **whitelist = read_txt(filename_whitelist, nrow);

    qsort(whitelist, nrow, sizeof(char *), compare);

    tree_whitelist = construct_tree(whitelist, 0, nrow - 1);

    // print_tree(tree_whitelist);

    int t1 = in(tree_whitelist, "TTTGGTTGTGACCAAG");

    printf("t1 = %d\n", t1);

    // print_tree(tree_whitelist);
    int *prod_id = malloc(sizeof(int));
    int *cons_id = malloc(NUM_CONSUMERS * sizeof(int));

    pthread_t producer_thread[NUM_PRODUCERS];
    pthread_t consumer_thread[NUM_CONSUMERS];

    for (int i = 0; i < NUM_PRODUCERS; i++)
    {
        prod_id[i] = i;
        pthread_create(&producer_thread[i], NULL, producer, &prod_id[i]);
    }

    for (int i = 0; i < NUM_CONSUMERS; i++)
    {
        cons_id[i] = i;
        pthread_create(&consumer_thread[i], NULL, consumer, &cons_id[i]);
    }

    // wait for the reader and processor threads to finish
    for (int i = 0; i < NUM_PRODUCERS; i++)
    {
        pthread_join(producer_thread[i], NULL);
    }

    for (int i = 0; i < NUM_CONSUMERS; i++)
    {

        pthread_join(consumer_thread[i], NULL);
    }

    free(cons_id);
    free(reader_buffer);

    free_tree_node(tree_whitelist);
    free(whitelist);

    for (int i = 0; i < 3; i++)
    {
        gzclose(file_in[i]);
        gzclose(file_out[i]);
    }

    return 0;
}