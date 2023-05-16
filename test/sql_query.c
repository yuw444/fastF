#include <stdio.h>
#include <stdlib.h>
#include <sqlite3.h>
#include <string.h>
#include <zlib.h>

int main(int argc, char *argv[])
{
    sqlite3 *db;
    char *zErrMsg = 0;
    int rc;
    char *sql;
    int cell_id;
    int gene_id;
    char UMI[20];

    rc = sqlite3_open("mydatabase.db", &db);

    if (rc)
    {
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return (1);
    }
    else
    {
        printf("Opened database successfully\n");
    }

    gzFile fp = gzopen("/scratch/u/yu89975/fastF/build/test/mtx.gz", "wb");

    if (fp == NULL)
    {
        fprintf(stderr, "Can't open file: %s\n", "/scratch/u/yu89975/fastF/build/test/mtx.gz");
        return (1);
    }

    sqlite3_stmt *stmt;
    // Prepare the query
    sqlite3_prepare_v2(db, "SELECT cell_id, gene_id, COUNT(DISTINCT UMI) as distinct_count FROM UMI GROUP BY cell_id, gene_id", -1, &stmt, NULL);
    // Execute the query
    while (sqlite3_step(stmt) == SQLITE_ROW)
    {
        // Get the values from the result row
        int id1 = sqlite3_column_int(stmt, 0);
        int id2 = sqlite3_column_int(stmt, 1);
        int distinct_count = sqlite3_column_int(stmt, 2);
        // printf("ID1: %d, ID2: %d, Distinct Count: %d\n", id1, id2, distinct_count);
        gzprintf(fp, "%d %d %d\n", id2, id1, distinct_count);
    }
    sqlite3_finalize(stmt);


    sqlite3_close(db);
    gzclose(fp);

    return (0);
}
