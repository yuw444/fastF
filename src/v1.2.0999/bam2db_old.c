#include "bam2db.h"

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

void bam2db(
    char *bam_file,
    char *db_file,
    char *cell_barcode_file,
    char *feature_name_file)
{
    sqlite3 *db;
    char *zErrMsg = 0;
    int rc;
    char *sql;
    sqlite3_stmt *stmt;

    // open database
    rc = sqlite3_open(db_file, &db);
    if (rc)
    {
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        exit(1);
    }
    else
    {
        fprintf(stderr, "Opened database successfully\n");
    }

    // open bam file
    samFile *bam = sam_open(bam_file, "r");
    if (bam == NULL)
    {
        fprintf(stderr, "Fail to open BAM file %s\n", bam_file);
        exit(1);
    }
    else
    {
        fprintf(stderr, "Opened BAM file %s successfully\n", bam_file);
    }

    // open cell barcode file
    gzFile cell_barcode = gzopen(cell_barcode_file, "r");
    if (cell_barcode == NULL)
    {
        fprintf(stderr, "Fail to open cell barcode file %s\n", cell_barcode_file);
        exit(1);
    }
    else
    {
        fprintf(stderr, "Opened cell barcode file %s successfully\n", cell_barcode_file);
    }

    // open feature name file
    gzFile feature_name = gzopen(feature_name_file, "r");
    if (feature_name == NULL)
    {
        fprintf(stderr, "Fail to open feature name file %s\n", feature_name_file);
        exit(1);
    }
    else
    {
        fprintf(stderr, "Opened feature name file %s successfully\n", feature_name_file);
    }

    // create cell table
    sql = "CREATE TABLE cell (cell_index INTEGER AUTO_INCREMENT PRIMARY KEY, cell_barcode TEXT);";
    rc = sqlite3_exec(db, sql, NULL, 0, &zErrMsg);

    if (rc != SQLITE_OK)
    {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
    }

    // create feature table
    sql = "CREATE TABLE feature (feature_index INTEGER AUTO_INCREMENT PRIMARY KEY, feature_name TEXT);";
    rc = sqlite3_exec(db, sql, NULL, 0, &zErrMsg);

    if (rc != SQLITE_OK)
    {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
    }

    // create umi table
    sql = "CREATE TABLE umi (cell_index INTEGER, feature_index INTEGER, encoded_umi TEXT);";
    rc = sqlite3_exec(db, sql, NULL, 0, &zErrMsg);

    if (rc != SQLITE_OK)
    {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
    }

    // create hash table for cell barcode
    hash_table *ht_cell = hash_table_create(1 << 20, hash, NULL);
    if (ht_cell == NULL)
    {
        fprintf(stderr, "ERROR: Cannot create hash table for cell barcode\n");
        exit(1);
    }

    // create hash table for gene
    hash_table *ht_feature = hash_table_create(1 << 20, hash, NULL);
    if (ht_feature == NULL)
    {
        fprintf(stderr, "ERROR: Cannot create hash table for feature\n");
        exit(1);
    }

    // read cell barcode file and insert into hash table and database
    char *cell_barcode_buffer = (char *)malloc(1024 * sizeof(char));
    int cell_index = 1;

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);
    sql = "INSERT INTO cell VALUES (?1, ?2);";
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    while (gzgets(cell_barcode, cell_barcode_buffer, 1024) != NULL)
    {
        cell_barcode_buffer[strcspn(cell_barcode_buffer, "\n\r\t")] = '\0';
        int *cell_index_ptr = (int *)calloc(1, sizeof(int));
        *cell_index_ptr = cell_index;
        if (hash_table_insert(ht_cell, cell_barcode_buffer, cell_index_ptr))
        {
            sqlite3_bind_int(stmt, 1, cell_index);
            sqlite3_bind_text(stmt, 2, cell_barcode_buffer, strlen(cell_barcode_buffer), SQLITE_STATIC);
            if (sqlite3_step(stmt) != SQLITE_DONE)
            {
                fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
            }
            sqlite3_reset(stmt);
            cell_index++;
        }
        else
        {
            printf("Warning: Duplicate cell barcodes were found in %s!\n", cell_barcode_file);
        }
    }
    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);
    sqlite3_finalize(stmt);
    gzclose(cell_barcode);
    free(cell_barcode_buffer);

    // read feature name file and insert into hash table and database

    char *feature_name_buffer = (char *)malloc(1024 * sizeof(char));
    int feature_index = 1;

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);
    sql = "INSERT INTO feature VALUES (?1, ?2);";
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    while (gzgets(feature_name, feature_name_buffer, 1024) != NULL)
    {
        feature_name_buffer[strcspn(feature_name_buffer, "\n\r\t")] = '\0';
        int *feature_index_ptr = (int *)calloc(1, sizeof(int));
        *feature_index_ptr = feature_index;
        if (hash_table_insert(ht_feature, feature_name_buffer, feature_index_ptr))
        {
            sqlite3_bind_int(stmt, 1, feature_index);
            sqlite3_bind_text(stmt, 2, feature_name_buffer, strlen(feature_name_buffer), SQLITE_STATIC);
            if (sqlite3_step(stmt) != SQLITE_DONE)
            {
                fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
            }
            sqlite3_reset(stmt);
            feature_index++;
        }
        else
        {
            printf("Warning: Duplicate feature names were found in %s!\n", feature_name_file);
        }

    }

    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);
    sqlite3_finalize(stmt);
    gzclose(feature_name);
    free(feature_name_buffer);

    /*********convert bam file to sqlite3 database***********/
    bam_hdr_t *bam_header = sam_hdr_read(bam);
    bam1_t *bam_record = bam_init1();
    size_t total_reads_counts = 0;
    size_t valid_reads_counts = 0;
    time_t t = time(NULL);
    struct tm *tm;
    char s[64];

    printf("Start to convert bam file to sqlite3 database...\n");
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);
    sql = "INSERT INTO umi VALUES (?1, ?2, ?3);";
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);

    while (sam_read1(bam, bam_header, bam_record) >= 0)
    {
        total_reads_counts++;

        if (total_reads_counts % 1000000 == 0)
        {
            tm = localtime(&t);
            strftime(s, sizeof(s), "%c", tm);
            printf("%s: Processed %zu reads...\n", s, total_reads_counts);
        }

        uint8_t *cell_barcode_ptr = bam_aux_get(bam_record, "CB");

        if (cell_barcode_ptr == NULL)
        {
            continue;
        }

        char *cell_barcode = bam_aux2Z(cell_barcode_ptr);
        void *temp1 = hash_table_lookup(ht_cell, cell_barcode);

        if (temp1 == NULL)
        {
            continue;
        }

        int cell_index = *(int *)temp1;
        uint8_t *xf_ptr = bam_aux_get(bam_record, "xf");
        int UMI_quality = bam_aux2i(xf_ptr);

        if (UMI_quality != 25)
        {
            continue;
        }

        uint8_t *gx_ptr = bam_aux_get(bam_record, "GX");
        char *feature_ID = bam_aux2Z(gx_ptr);
        void *temp2 = hash_table_lookup(ht_feature, feature_ID);

        if (temp2 == NULL)
        {
            continue;
        }

        int feature_index = *(int *)temp2;
        uint8_t *ub_ptr = bam_aux_get(bam_record, "UB");
        char *UMI = bam_aux2Z(ub_ptr);
        // printf("%s\n", UMI);
        uint8_t *encoded_UMI = encode_DNA(UMI);
        // printf("%s\n", encoded_UMI);
        // char *decoded_UMI = decode_DNA(encoded_UMI, 10);
        // printf("%s\n", decoded_UMI);

        sqlite3_bind_int(stmt, 1, cell_index);
        sqlite3_bind_int(stmt, 2, feature_index);
        sqlite3_bind_text(stmt, 3, encoded_UMI, strlen(encoded_UMI), SQLITE_STATIC);

        if (sqlite3_step(stmt) != SQLITE_DONE)
        {
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
        }

        sqlite3_reset(stmt);
        valid_reads_counts++;

        free(encoded_UMI);
    }

    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);
    sqlite3_finalize(stmt);

    printf("Total reads: %zu\n", total_reads_counts);
    printf("Valid reads: %zu\n", valid_reads_counts);

    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    hts_close(bam);

    // create index for umi table for speeding up the query
    printf("Create index for umi table...\n");
    sql = "CREATE INDEX umi_index ON umi (cell_index, feature_index);";
    sqlite3_exec(db, sql, NULL, NULL, &zErrMsg);

    sqlite3_close(db);

    hash_table_destroy(ht_cell);
    hash_table_destroy(ht_feature);
}
