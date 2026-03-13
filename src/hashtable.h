#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

struct _entry {
    char *key;
    void *obj;
    struct _entry *next;
};

typedef uint64_t hashfunction(const char*, size_t);
typedef void cleanupfunction(void*);
typedef struct _entry entry;

struct _hash_table {
    uint32_t size;
    hashfunction *hash;
    cleanupfunction *cleanup;
    entry **elements;
};

typedef struct _hash_table hash_table;

hash_table *hash_table_create(uint32_t size, hashfunction *hf, cleanupfunction *cf);
void hash_table_destroy(hash_table *ht);
void hash_table_print(hash_table *ht);
// return true if the object was inserted, false if it already existed
bool hash_table_insert(hash_table *ht, const char *key, void *obj);
// return the object that was found
void *hash_table_lookup(hash_table *ht, const char *key);
// return the object that was deleted from the table, or NULL if it didn't exist
// user is responsible for freeing the object
void *hash_table_delete(hash_table *ht, const char *key);


#endif