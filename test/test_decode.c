#include <criterion/criterion.h>
#include "../src/v1.2.0999/bam2db.h"

Test (code, code_DNA)
{
    char *DNA_Seq1 = "ATCGATCGATGCTACC";
    char *DNA_Seq2 = "ATCGATCGATGCTACT";

    uint8_t *encoded_DNA1 = encode_DNA(DNA_Seq1);
    uint8_t *encoded_DNA2 = encode_DNA(DNA_Seq2);  

    if(strcmp(encoded_DNA1, encoded_DNA2) == 0) {
        printf("encoded_DNA1 == encoded_DNA2\n");
    } else {
        printf("encoded_DNA1 != encoded_DNA2\n");
    }

    printf("encoded_DNA1: %s\n", encoded_DNA1);
    printf("encoded_DNA2: %s\n", encoded_DNA2);

    char *decoded_DNA1 = decode_DNA(encoded_DNA1, strlen(DNA_Seq1));
    char *decoded_DNA2 = decode_DNA(encoded_DNA2, strlen(DNA_Seq2));

    printf("decoded_DNA1: %s\n", decoded_DNA1);
    printf("decoded_DNA2: %s\n", decoded_DNA2);

    cr_assert_str_eq(DNA_Seq1, decoded_DNA1);

    free(encoded_DNA1);
    free(encoded_DNA2);
    free(decoded_DNA1);
    free(decoded_DNA2);
}

Test (code, db)
{
    char *bam_file = "/scratch/u/yu89975/fastF_test/rawfastq/sub/AI/outs/possorted_genome_bam_filtered.bam";
    char *db_file = "/scratch/u/yu89975/fastF/data/test_pure.db";
    char *barcode_file = "/scratch/u/yu89975/fastF_test/rawfastq/sub/AI/outs/filtered_feature_bc_matrix/barcodes.tsv.gz";
    char *feature_file = "/scratch/u/yu89975/fastF_test/rawfastq/sub/AI/outs/filtered_feature_bc_matrix/features.tsv.gz";

    bam2db(bam_file, db_file, barcode_file, feature_file);
}