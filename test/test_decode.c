#include <criterion/criterion.h>
#include "../src/bam2db_ds.h"

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

Test(code, end)
{
    char temp[] = "ATCGATCGATGCTACC";
    temp[4] = '\0';

    printf("%s\n", temp);
    printf("%s\n", temp+5);

    char *temp2 = (char *)malloc(18 * sizeof(char));
    memcpy(temp2, temp, 17);

    printf("%s\n", temp2);
    printf("%s\n", temp2+5);

    if (strncmp(temp, temp2, 18) == 0) {
        printf("temp == temp2\n");
    } else {
        printf("temp != temp2\n");
    }

    free(temp2);
}