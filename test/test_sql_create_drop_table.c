#include <stdio.h>
#include <stdlib.h>
#include <sqlite3.h>

size_t nrow_sql_table (
    sqlite3 *db_handle,
    const char *table_name)
{
    char *err_msg = NULL;
    int rc;
    sqlite3_stmt *stmt;

    char *sql = (char *)malloc(1024 * sizeof(char));
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

int main() {
    sqlite3* db;
    char* errorMessage = NULL;
    const char* databasePath = "/scratch/u/yu89975/fastF/test/mydatabase.db";

    // Open the database
    int result = sqlite3_open(databasePath, &db);
    if (result != SQLITE_OK) {
        fprintf(stderr, "Error opening the database: %s\n", sqlite3_errmsg(db));
        return 1;
    }

    // SQL statement to create a table
    const char* createTableSQL = "CREATE TABLE IF NOT EXISTS UMIdss AS "
                                 "SELECT * " 
                                 "FROM UMI "
                                 "LIMIT 10000;";

    // Execute the create table statement
    result = sqlite3_exec(db, createTableSQL, NULL, NULL, &errorMessage);
    if (result != SQLITE_OK) {
        fprintf(stderr, "Error creating table: %s\n", errorMessage);
        sqlite3_free(errorMessage);
        sqlite3_close(db);
        return 1;
    }

    printf("Table created successfully.\n");

    // SQL statement to drop the table
    const char* dropTableSQL = "DROP TABLE IF EXISTS UMIds;";

    // Execute the drop table statement
    result = sqlite3_exec(db, dropTableSQL, NULL, NULL, &errorMessage);
    if (result != SQLITE_OK) {
        fprintf(stderr, "Error dropping table: %s\n", errorMessage);
        sqlite3_free(errorMessage);
        sqlite3_close(db);
        return 1;
    }

    printf("Table dropped successfully.\n");

    size_t nrows = nrow_sql_table(db, "UMIdss");

    printf("Number of rows: %zu\n", nrows);

    // Close the database connection
    sqlite3_close(db);
    return 0;

}
