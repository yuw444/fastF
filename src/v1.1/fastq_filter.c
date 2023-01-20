// Author : Yu Wang
// Date : 01/13/2023
/*

This program is used to filter the file files based on the whitelist.
It is designed to be run with multiple threads.

*/

#include "fastq_filter.h"

long int num_waits = 0;
int flag = 0; // flag to indicate the end of the program

pthread_mutex_t queue_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t check_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t write_mutex = PTHREAD_MUTEX_INITIALIZER;

pthread_cond_t not_empty = PTHREAD_COND_INITIALIZER;
pthread_cond_t not_full = PTHREAD_COND_INITIALIZER;

queue *init_queue()
{
    queue *q = (queue *)malloc(sizeof(struct queue));
    q->head = 0;
    q->tail = 0;
    q->size = 0;
    return q;
}

bool is_empty(struct queue *q)
{
    return q->size == 0;
}

bool is_full(struct queue *q)
{
    return q->size == MAX_QUEUE_SIZE;
}

void free_fastq(fastq *block)
{
    free(block->id);
    free(block->seq);
    free(block->qual);
    free(block);
}
void free_comb_fastq(comb_fastq *comb)
{
    free_fastq(comb->I1);
    free_fastq(comb->R1);
    free_fastq(comb->R2);

    free(comb);
}

void enqueue(struct queue *q, comb_fastq *value)
{
    pthread_mutex_lock(&queue_lock);
    while (is_full(q))
    {
        printf("queue is full, waiting deque!\n");
        pthread_cond_wait(&not_full, &queue_lock);
    }
    value->random_number = (double)rand() / (double)RAND_MAX;
    q->data[q->tail] = value;
    q->tail = (q->tail + 1) % MAX_QUEUE_SIZE;
    q->size++;
    pthread_cond_signal(&not_empty);
    pthread_mutex_unlock(&queue_lock);
}

comb_fastq *dequeue(struct queue *q)
{
    pthread_mutex_lock(&queue_lock);
    while (is_empty(q))
    {
        printf("queue is empty, waiting enque! %ld \n", num_waits++);
        pthread_cond_wait(&not_empty, &queue_lock);
    }
    comb_fastq *value = q->data[q->head];
    q->head = (q->head + 1) % MAX_QUEUE_SIZE;
    q->size--;
    pthread_cond_signal(&not_full);
    pthread_mutex_unlock(&queue_lock);

    return value;
}

fastq *get_fastq(gzFile file)
{
    fastq *out = (fastq *)malloc(sizeof(fastq));
    out->id = (char *)malloc(MAX_LINE_LENGTH);
    out->seq = (char *)malloc(MAX_LINE_LENGTH);
    out->qual = (char *)malloc(MAX_LINE_LENGTH);

    char temp[MAX_LINE_LENGTH];
    if (gzgets(file, out->id, MAX_LINE_LENGTH) != NULL)
    {
        gzgets(file, out->seq, MAX_LINE_LENGTH);
        gzgets(file, temp, MAX_LINE_LENGTH);
        gzgets(file, out->qual, MAX_LINE_LENGTH);
    }
    else
    {
        flag = 1;
        free_fastq(out);
        return NULL;
    }

    return out;
}

comb_fastq *get_comb_fastq(gzFile fastq[3])
{
    comb_fastq *out = (comb_fastq *)malloc(sizeof(comb_fastq));
    out->I1 = get_fastq(fastq[0]);
    out->R1 = get_fastq(fastq[1]);
    out->R2 = get_fastq(fastq[2]);

    if (out->I1->id == NULL || out->R1->id == NULL || out->R2->id == NULL)
    {
        free_comb_fastq(out);
        return NULL;
    }
    return out;
}


// sorting function for whitelist

int compare(const void *a, const void *b)
{
    return strcmp(*(const char **)a, *(const char **)b);
}

// construct a binary searching tree for the whitelist
node *construct_tree(char **sorted_string, int start, int end)
{
    if (start > end)
    {
        return NULL;
    }

    int mid = (start + end) / 2;
    node *root = (node *)malloc(sizeof(node));
    root->data = sorted_string[mid];
    root->left = construct_tree(sorted_string, start, mid - 1);
    root->right = construct_tree(sorted_string, mid + 1, end);

    return root;
}

// print the binary searching tree of whitelist
void print_tree(node *root)
{
    if (root == NULL)
    {
        return;
    }
    printf("%s", root->data);
    print_tree(root->left);
    print_tree(root->right);
}

// free binary tree node

void free_tree_node(node *root)
{
    if (root == NULL)
    {
        return;
    }
    free(root->data);
    free_tree_node(root->left);
    free_tree_node(root->right);
    free(root);
}

// get row number of text file

int get_row(char *file_name)
{
    int i = 0;
    char buffer[MAX_LINE_LENGTH];
    FILE *stream = fopen(file_name, "r");
    if (stream == NULL)
    {
        printf("Error opening %s!\n", file_name);
        exit(1);
    }

    while (fgets(buffer, 1024, stream) != NULL)
    {
        i++;
        // printf("%s", buffer);
    }
    fclose(stream);

    return i;
}


// read whitelist from file

char **read_txt(char *file_name, size_t nrows)
{
    FILE *stream = fopen(file_name, "r");

    if (stream == NULL)
    {
        printf("Error opening %s!\n", file_name);
        exit(1);
    }

    char **whitelist = (char **)malloc(nrows * sizeof(char *));
    int i = 0;

    for (int i = 0; i < nrows; i++)
    {
        whitelist[i] = (char *)malloc((LEN_CELLBARCODE + 2) * sizeof(char));
        fgets(whitelist[i], LEN_CELLBARCODE + 2, stream);
        // printf("%s %ld", whitelist[i], strlen(whitelist[i]));
    }

    fclose(stream);

    return whitelist;
}

bool in(node *root, char *element)
{
    if (root == NULL)
    {
        return 0;
    }

    if (strncmp(root->data, element, LEN_CELLBARCODE) == 0)
    {
        return 1;
    }
    else if (strncmp(root->data, element, LEN_CELLBARCODE) > 0)
    {
        return in(root->left, element);
    }
    else
    {
        return in(root->right, element);
    }
}

// get substring of a string

char *substring(char *string, int position, int length)
{
    char *pointer =(char *) malloc(length + 1);

    if (pointer == NULL)
    {
        printf("Unable to allocate memory.\n");
        exit(1);
    }

    strncpy(pointer, string + position, length);

    *(pointer + length) = '\0';

    return pointer;
}

// combine struct fastq to one whole string

char *combine_string(fastq *block)
{
    char *buf = (char *)malloc(MAX_LINE_LENGTH * 3 * sizeof(char));
    snprintf(buf, MAX_LINE_LENGTH * 3, "%s%s+\n%s", block->id, block->seq, block->qual);
    return buf;
}

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

            // printf("random number: %f \n", comb->random_number);

            if (comb->random_number <= rate_threshold && in(tree_whitelist, cell_barcode))
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

// int main()
// {
//     double startTime = clock();

//     unsigned int seed = 926;
//     srand(seed);

//     file_in[0] = gzopen("/home/rstudio/Frag/data/AI_S1_L001_I1_001.fastq.gz", "rb");
//     file_in[1] = gzopen("/home/rstudio/Frag/data/AI_S1_L001_R1_001.fastq.gz", "rb");
//     file_in[2] = gzopen("/home/rstudio/Frag/data/AI_S1_L001_R1_001.fastq.gz", "rb");

//     file_out[0] = gzopen("/home/rstudio/Frag/data/I1_thread.fastq.gz", "wb");
//     file_out[1] = gzopen("/home/rstudio/Frag/data/R1_thread.fastq.gz", "wb");
//     file_out[2] = gzopen("/home/rstudio/Frag/data/R2_thread.fastq.gz", "wb");

//     reader_buffer = init_queue();

//     rate_threshold = 0.8;

//     char *filename_whitelist = "../data/whitelist.txt";

//     printf("Reading whitelist...\n");
//     int nrow = get_row(filename_whitelist);

//     printf("nrow = %d\n", nrow);

//     char **whitelist = read_txt(filename_whitelist, nrow);

//     qsort(whitelist, nrow, sizeof(char *), compare);

//     tree_whitelist = construct_tree(whitelist, 0, nrow - 1);

//     // print_tree(tree_whitelist);

//     int t1 = in(tree_whitelist, "TTTGGTTGTGACCAAG");

//     printf("t1 = %d\n", t1);

//     // print_tree(tree_whitelist);
//     int *prod_id = malloc(NUM_PRODUCERS * sizeof(int));
//     int *cons_id = malloc(NUM_CONSUMERS * sizeof(int));

//     pthread_t producer_thread[NUM_PRODUCERS];
//     pthread_t consumer_thread[NUM_CONSUMERS];

//     for (int i = 0; i < NUM_PRODUCERS; i++)
//     {
//         prod_id[i] = i;
//         pthread_create(&producer_thread[i], NULL, producer, &prod_id[i]);
//     }

//     sleep(1);

//     for (int i = 0; i < NUM_CONSUMERS; i++)
//     {
//         cons_id[i] = i;
//         pthread_create(&consumer_thread[i], NULL, consumer, &cons_id[i]);
//     }

//     // wait for the reader and processor threads to finish
//     for (int i = 0; i < NUM_PRODUCERS; i++)
//     {
//         pthread_join(producer_thread[i], NULL);
//     }

//     for (int i = 0; i < NUM_CONSUMERS; i++)
//     {

//         pthread_join(consumer_thread[i], NULL);
//     }

//     free(cons_id);
//     free(reader_buffer);

//     free_tree_node(tree_whitelist);
//     free(whitelist);

//     for (int i = 0; i < 3; i++)
//     {
//         gzclose(file_in[i]);
//         gzclose(file_out[i]);
//     }

//     double stopTime = clock();

//     double secsElapsed = (stopTime - startTime) / CLOCKS_PER_SEC;

//     printf("Elapsed: %f seconds\n", secsElapsed);
//     return 0;
// }

// // Elapsed: 403.236648 seconds