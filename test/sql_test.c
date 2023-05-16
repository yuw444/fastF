#include <stdio.h>
#include <stdlib.h>
#include <sqlite3.h>
#include <string.h>
#include <zlib.h>

int main(int argc, char* argv[]){
    sqlite3 *db;
    char *zErrMsg = 0;
    int rc;
    char *sql;
    int cell_id;
    int gene_id;
    char UMI[20];

    rc = sqlite3_open("mydatabase.db", &db);

    if(rc){
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return(1);
    }else{
        printf("Opened database successfully\n");
    }

    gzFile fp = gzopen("/scratch/u/yu89975/fastF/build/test/rst.gz", "r");

    if(fp == NULL){
        fprintf(stderr, "Can't open file: %s\n", "/scratch/u/yu89975/fastF/build/test/test.rst.csv.gz");
        return(1);
    }

    // check if table exists, if so delete it
    sql = sqlite3_mprintf("DROP TABLE IF EXISTS UMI");
    rc = sqlite3_exec(db, sql, NULL, NULL, NULL);

    if(rc != SQLITE_OK){
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        // sqlite3_free(zErrMsg);
    }

    sql = "CREATE TABLE UMI("  \
          "cell_id  INT     NOT NULL," \
          "gene_id  INT     NOT NULL," \
          "UMI      TEXT    NOT NULL);";

    /* Execute SQL statement */
    rc = sqlite3_exec(db, sql, 0, 0, &zErrMsg);

    if(rc != SQLITE_OK){
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        // sqlite3_free(zErrMsg);
    }

    char buffer[1024];

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);
    char *sql_insert = "INSERT INTO UMI (cell_id, gene_id, UMI) VALUES (?1, ?2, ?3)";
    sqlite3_stmt *stmt;
    sqlite3_prepare_v2(db, sql_insert, strlen(sql_insert), &stmt, NULL);

    while(gzgets(fp, buffer, 1024) != Z_NULL){
        /* Create SQL statement */
        sscanf(buffer, "%d,%d,%s", &cell_id, &gene_id, UMI);

        sqlite3_bind_int(stmt, 1, cell_id);
        sqlite3_bind_int(stmt, 2, gene_id);
        sqlite3_bind_text(stmt, 3, UMI, strlen(UMI), SQLITE_STATIC);

        if (sqlite3_step(stmt) != SQLITE_DONE) {
            fprintf(stderr, "SQL error: %s\n", sqlite3_errmsg(db));
            // sqlite3_free(zErrMsg);
        }

        sqlite3_reset(stmt);

    }

    sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &zErrMsg);
    sqlite3_finalize(stmt);

    // sqlite3_free(sql);
    gzclose(fp);
    sqlite3_close(db);
    return(0);
}
