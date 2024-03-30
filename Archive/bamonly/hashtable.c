#include "hashtable.h"

struct _entry {
    char *key;
    void *obj;
    struct _entry *next;
};

struct _hash_table {
    uint32_t size;
    hashfunction *hash;
    cleanupfunction *cleanup;
    entry **elements;
};


static size_t hash_table_index(hash_table *ht, const char *key)
{
    size_t result = ht->hash(key, strlen(key)) % ht->size;
    return result;
}

// pass in NULL cf for default free behavior
hash_table *hash_table_create(uint32_t size, hashfunction *hf, cleanupfunction *cf)
{
    hash_table *ht = (hash_table *)malloc(sizeof(hash_table));
    ht->size = size;
    ht->hash = hf;
    if (cf) 
    {
        ht->cleanup = cf;
    }
    else
    {
        ht->cleanup = free;
    }
    // note that calloc zeros out the memory
    ht->elements = (entry **)calloc(size, sizeof(entry *));
    return ht;
}


void hash_table_destroy(hash_table *ht)
{

    for (uint32_t i = 0; i < ht->size; i++) {
        entry *e = ht->elements[i];
        while (e != NULL) {
            entry *next = e->next;
            free(e->key);
            ht->cleanup(e->obj);
            free(e);
            e = next;
        }
    }

    free(ht->elements);
    free(ht);
}


void hash_table_print(hash_table *ht)
{
    printf("Start printing hash table:\n");
    for (int i = 0; i < ht->size; i++) {
        entry *e = ht->elements[i];
        if (e == NULL) {
            // printf("\t%i\t ---\n", i);
        } else {
            printf("\t%i\t", i);
            while (e != NULL) {
                printf(" \"%s\"(%zu) -  ", e->key, *(size_t *)e->obj);
                e = e->next;
            }
            printf("\n");
        }
    }
    printf("End table\n");
}

// if the key already exists, return false and free the key and obj is user's responsibility
// if the key does not exist, return true and free is taken cared by the hash table destory function
bool hash_table_insert(hash_table *ht, const char *key, void *obj)
{
    if (ht == NULL || key == NULL || obj == NULL) return false;
    size_t index = hash_table_index(ht, key);

    // check if the key already exists, if so, return false
    if (hash_table_lookup(ht, key) != NULL)
    {
        // key already exists, and free the key and obj is user's responsibility
        return false;
    } 

    // create a new entry
    entry *e = (entry *)malloc(sizeof(entry));
    e->obj = obj;
    e->key = strdup(key);

    // insert the entry into the current head of index in the hash table
    // ht->elements[index] is the head of the linked list
    e->next = ht->elements[index];
    // update the head of the linked list
    ht->elements[index] = e;

    return true;

}


void *hash_table_lookup(hash_table *ht, const char *key)
{
    if(key == NULL || ht == NULL) return false;
    size_t index = hash_table_index(ht, key);
    // printf("index: %zu\n", index);

    entry *tmp = ht->elements[index];
    while (tmp != NULL && strcmp(tmp->key , key) != 0)
    {
        tmp = tmp->next;
    } 

    // if tmp is NULL, then the key does not exist
    if (tmp == NULL) return NULL;
        
    return tmp->obj;
    
}


void *hash_table_delete(hash_table *ht, const char *key)
{
    if (ht == NULL || key == NULL) return false;

    size_t index = hash_table_index(ht, key);

    entry *tmp = ht->elements[index];
    entry *prev = NULL;
    while (tmp != NULL && strcmp(tmp->key , key) != 0)
    {
        prev = tmp;
        tmp = tmp->next;
    }

    // if tmp is NULL, then the key does not exist
    if (tmp == NULL) return NULL;

    // if prev is NULL, then the key is the head of the linked list
    if (prev == NULL) {
        ht->elements[index] = tmp->next;
    } else {
        // otherwise, update the previous entry to point to the next entry; namely, skip over tmp
        prev->next = tmp->next;
    }

    void *result = tmp->obj;
    free(tmp);
    return result;
}