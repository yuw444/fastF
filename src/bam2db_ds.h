#ifndef BAM2DB_H
#define BAM2DB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <htslib/sam.h>
#include <zlib.h>
#include <sqlite3.h>
#include <sys/stat.h>
#include "hashtable.h"
#include "mt19937ar.h"
#include "utils.h"

#define MAX_LINE_LENGTH 1024
#define BASE_BITS 2                      // the number of bits to encode a base
#define BYTE_SIZE 8                      // the number of bits in a byte
#define BYTE_MASK 0xFF                   // the mask of a byte
#define BASE_MASK ((1 << BASE_BITS) - 1) // the mask of a base

extern int _umi_copies_flag;

/**
 * @brief encode a base to uint8_t
 * 
 * @param base one of A, C, G, T
 * @return uint8_t 
 */
uint8_t encode_base(char base);

/**
 * @brief encode a DNA sequence to uint8_t pointer, the DNA sequence must be '\0' terminated
 * 
 * @param DNA_seq DNA sequence
 * @return uint8_t* 
 */
uint8_t *encode_DNA(const char *DNA_seq);

/**
 * @brief decode to char* with length n from unit8_t
 * 
 * @param encoded_DNA  the encoded DNA sequence
 * @param n the length of the encoded DNA sequence
 * @return char* 
 */
char *decode_DNA(uint8_t *encoded_DNA, size_t n);

/**
 * @brief convert bam file to sql db file
 * 
 * @param bam_file path to bam file
 * @param db_file path to db file
 * @param cell_barcode_file path to cell barcode file
 * @param feature_name_file path to feature name file
 * @param rate_cell the rate of cell barcode
 * @param rate_depth the rate of depth
 * @param seed the seed of random
 * @return 0 if success, 1 if fail
 */
int bam2db(
    char *bam_file,
    char *db_file,
    char *path_out,
    char *barcodes_file,
    char *features_file,
    float rate_cell,
    float rate_depth,
    unsigned int seed);



/**
 * @brief convert sql table to gz file
 * 
 * @param db_handle sqlite3 handle
 * @param table_name table name
 * @param gz_file_ptr gz file pointer
 * @param header whether to write header
 * @param delim delimiter
 * 
 * @return 0 if success, 1 if fail
 */

int table2gz(
    sqlite3 *db_handle,
    const char *table_name, 
    gzFile gz_file_ptr,
    unsigned int header, 
    const char *delim);

/**
 * @brief  get the number of rows in a sql table
 * 
 * @param db_handle sqlite3 handle
 * @param table_name table name
 * 
 * @return the number of rows in a sql table
 */

size_t nrow_sql_table (
    sqlite3 *db_handle,
    const char *table_name);
#endif
