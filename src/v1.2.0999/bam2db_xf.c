#include "bam2db_ds.h"

uint8_t encode_base(char base)
{
    switch (base)
    {
    case 'A':
        return 0b00;
    case 'C':
        return 0b01;
    case 'G':
        return 0b10;
    case 'T':
        return 0b11;
    default:
        return BYTE_MASK; // unknown base
    }
}

uint8_t *encode_DNA(const char *DNA_seq)
{
    uint16_t len_DNA = strlen(DNA_seq);
    // plus 1 for the end of string if like to print it out
    uint8_t *encoded_DNA = (uint8_t *)calloc((len_DNA + BYTE_SIZE / BASE_BITS - 1) / 4 + 1, sizeof(uint8_t));
    if (encoded_DNA == NULL)
    {
        return NULL;
    }
    uint8_t *p = encoded_DNA;
    int bit_index = BYTE_SIZE - BASE_BITS;
    for (size_t i = 0; i < len_DNA; i++)
    {
        uint8_t base = encode_base(DNA_seq[i]);
        if (base == BYTE_MASK)
        {
            free(encoded_DNA);
            return NULL; // unknown base
        }
        *p |= (base & BASE_MASK) << bit_index;
        bit_index -= BASE_BITS;
        if (bit_index < 0)
        {
            bit_index = BYTE_SIZE - BASE_BITS;
            p++;
        }
    }
    return encoded_DNA;
}

char *decode_DNA(uint8_t *encoded_DNA, size_t n)
{
    char *decoded_DNA = (char *)calloc(n + 1, sizeof(char));
    decoded_DNA[n] = '\0';
    if (decoded_DNA == NULL)
    {
        return NULL;
    }
    char *p = decoded_DNA;
    int bit_index = BYTE_SIZE - BASE_BITS;
    for (size_t i = 0; i < n; i++)
    {
        uint8_t base = (encoded_DNA[i / 4] >> bit_index) & BASE_MASK;
        switch (base)
        {
        case 0b00:
            *p = 'A';
            break;
        case 0b01:
            *p = 'C';
            break;
        case 0b10:
            *p = 'G';
            break;
        case 0b11:
            *p = 'T';
            break;
        default:
            free(decoded_DNA);
            return NULL; // unknown base
        }
        bit_index -= BASE_BITS;
        if (bit_index < 0)
        {
            bit_index = BYTE_SIZE - BASE_BITS;
        }
        p++;
    }
    return decoded_DNA;
}

// hash function
uint64_t hash(const char *str, size_t len)
{
    uint64_t hash = 5381;

    for (int i = 0; i < len; i++)
        hash = ((hash << 5) + hash) + str[i]; /* hash * 33 + c */

    return hash;
}


int table2gz(
    sqlite3 *db_handle,
    const char *table_name, 
    gzFile gz_file_ptr,
    unsigned int header, 
    const char *delim)
{
    char *sql = malloc(1024 * sizeof(char));

    // prepare the sql statement
    sprintf(sql, "SELECT * FROM %s;", table_name);
    sqlite3_stmt *stmt;
    if (sqlite3_prepare_v2(db_handle, sql, -1, &stmt, NULL) != SQLITE_OK)
    {
        fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db_handle));
        gzclose(gz_file_ptr);
        return EXIT_FAILURE;
    }

    // get the number of columns
    int column_count = sqlite3_column_count(stmt);

    if (header == 1) 
    {
        for (int i = 0; i < column_count; i++)
        {
            gzprintf(gz_file_ptr, "%s", sqlite3_column_name(stmt, i));
            if (i < column_count - 1)
            {
                gzprintf(gz_file_ptr, "%s", delim);
            }
        }

        gzprintf(gz_file_ptr, "\n");
    }

    // loop through the rows
    while (sqlite3_step(stmt) == SQLITE_ROW)
    {
        for (int i = 0; i < column_count; i++)
        {
            switch (sqlite3_column_type(stmt, i))
            {
                case SQLITE_INTEGER:
                    gzprintf(gz_file_ptr, "%d", sqlite3_column_int(stmt, i));
                    break;
                case SQLITE_FLOAT:
                    gzprintf(gz_file_ptr, "%.3f", sqlite3_column_double(stmt, i));
                    break;
                case SQLITE_TEXT:
                    gzprintf(gz_file_ptr, "%s", sqlite3_column_text(stmt, i));
                    break;
                default:
                    gzprintf(gz_file_ptr, "%s", "NULL");
                    break;
            }
            if (i < column_count - 1)
            {
                gzprintf(gz_file_ptr, "%s", delim);
            }
        }
        gzprintf(gz_file_ptr, "\n");
    }

    sqlite3_finalize(stmt);
    // gzclose(gz_file_ptr);
    free(sql);
    return EXIT_SUCCESS;
}

// get the number of rows in a table of sqlite3 database
size_t nrow_sql_table (
    sqlite3 *db_handle,
    const char *table_name)
{
    char *err_msg = NULL;
    int rc;
    sqlite3_stmt *stmt;

    char *sql = malloc(1024 * sizeof(char));
    sprintf(sql, "SELECT COUNT(*) FROM %s;", table_name);

    rc = sqlite3_prepare_v2(db_handle, sql, -1, &stmt, NULL);

    if (rc != SQLITE_OK)
    {
        fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db_handle));
        sqlite3_close(db_handle);
        return 0;
    }

    rc = sqlite3_step(stmt);

    if (rc != SQLITE_ROW)
    {
        fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db_handle));
        sqlite3_close(db_handle);
        return 0;
    }

    size_t nrow = sqlite3_column_int(stmt, 0);

    sqlite3_finalize(stmt);

    free(sql);

    return nrow;

}

