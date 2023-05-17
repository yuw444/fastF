#include <stdio.h>
#include <stdlib.h>
#include <sqlite3.h>
#include <zlib.h>

int writeTableToGzFile(sqlite3* db, const char* tableName, const char* gzFilePath) {
    gzFile gzFilePtr = gzopen(gzFilePath, "wb");
    if (gzFilePtr == NULL) {
        fprintf(stderr, "Error opening the Gzipped file for writing\n");
        return 1;
    }

    // Prepare the SQL statement
    char query[256];
    snprintf(query, sizeof(query), "SELECT gene_id, cell_id, COUNT(DISTINCT UMI) AS expression_level FROM %s GROUP BY cell_id, gene_id;", tableName);

    // Create a prepared statement from the query
    sqlite3_stmt* stmt;
    if (sqlite3_prepare_v2(db, query, -1, &stmt, NULL) != SQLITE_OK) {
        fprintf(stderr, "Error preparing SQLite statement: %s\n", sqlite3_errmsg(db));
        gzclose(gzFilePtr);
        return 1;
    }

    // Get the number of columns in the result set
    int columnCount = sqlite3_column_count(stmt);

    // Write column names to the Gzipped file
    for (int i = 0; i < columnCount; i++) {
        const char* columnName = sqlite3_column_name(stmt, i);
        gzprintf(gzFilePtr, "%s", columnName);
        if (i < columnCount - 1)
            gzputc(gzFilePtr, '\t');
    }
    gzputc(gzFilePtr, '\n');

    // Fetch and write table data row by row
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        for (int i = 0; i < columnCount; i++) {
            if (sqlite3_column_type(stmt, i) == SQLITE_INTEGER) {
                int columnValue = sqlite3_column_int(stmt, i);
                gzprintf(gzFilePtr, "%d", columnValue);
            } else {
                const unsigned char* columnValue = sqlite3_column_text(stmt, i);
                gzprintf(gzFilePtr, "%s", columnValue);
            }
            if (i < columnCount - 1)
                gzputc(gzFilePtr, '\t');
        }
        gzputc(gzFilePtr, '\n');
    }
    sqlite3_finalize(stmt);
    gzclose(gzFilePtr);
    return 0;
}

int main() {
    sqlite3* db;
    const char* databasePath = "/scratch/u/yu89975/fastF/test/mydatabase.db";
    const char* tableName = "UMIds";
    const char* gzFilePath = "output.gz";

    int result = sqlite3_open(databasePath, &db);
    if (result != SQLITE_OK) {
        fprintf(stderr, "Error opening the database: %s\n", sqlite3_errmsg(db));
        return 1;
    }

    if (writeTableToGzFile(db, tableName, gzFilePath) == 0)
        printf("Table '%s' has been written to the Gzipped file '%s'.\n", tableName, gzFilePath);
    else
        printf("Error writing the table '%s' to the Gzipped file.\n", tableName);

    sqlite3_close(db);
    return 0;
}
