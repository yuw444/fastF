
#include "hashtable.h"
#define MAX_LINE_LENGTH 1024

uint64_t hash(const char *str, size_t len)
{
    uint64_t hash = 5381;
    
    for (int i=0; i < len; i++)
        hash = ((hash << 5) + hash) + str[i]; /* hash * 33 + c */

    return hash;
}

void mycleanup(void *obj)
{
    free(obj);
}

int main(int argc, char **argv)
{
    if(argc != 3) 
    {
        printf("Usage: %s <wordlist> <number guesses>\n", argv[0]);
        return EXIT_FAILURE;
    }

    char *wordlist = argv[1];
    uint32_t num_guesses = atol(argv[2]);

    const int tablesize = (1<<20);
    hash_table *ht = hash_table_create(tablesize, hash, NULL);

    FILE *fp = fopen(wordlist, "r");
    if(fp == NULL)
    {
        printf("Error opening file %s\n", wordlist);
        return EXIT_FAILURE;
    }

    uint32_t num_words = 0;
    char buffer[MAX_LINE_LENGTH];

    while (!feof(fp) && fgets(buffer, MAX_LINE_LENGTH, fp) != NULL)
    {
        buffer[strcspn(buffer, "\n\r")] = '\0';
        char *newentry = (char *)malloc(strlen(buffer) + 1);
        strcpy(newentry, buffer);
        if(!hash_table_insert(ht, newentry, newentry))
        {
            free(newentry);
        }
        num_words++;
    }

    fclose(fp);
    printf("Read %d words from %s\n", num_words, wordlist);

    // hash_table_print(ht);
    hash_table_destroy(ht);


}