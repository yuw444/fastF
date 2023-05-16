#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define BASE_BITS 2 // the number of bits to encode a base
#define BYTE_SIZE 8 // the number of bits in a byte
#define BYTE_MASK 0xFF // the mask of a byte
#define BASE_MASK ((1 << BASE_BITS) - 1) // the mask of a base

uint8_t encode_base(char base) {
    switch (base) {
        case 'A': return 0b00;
        case 'C': return 0b01;
        case 'G': return 0b10;
        case 'T': return 0b11;
        default: return BYTE_MASK; // unknown base
    }
}

uint8_t* encode_UMI(const char* UMI) {

    if(strlen(UMI) != 10) {
        printf("Error: UMI length is not 10\n");
        exit(1);
    }

    uint8_t* encoded_UMI = (uint8_t*)calloc(4, sizeof(uint8_t));
    if (encoded_UMI == NULL) {
        return NULL;
    }
    uint8_t* p = encoded_UMI;
    int bit_index = BYTE_SIZE - BASE_BITS;
    for (size_t i = 0; i < 10; i++) {
        uint8_t base = encode_base(UMI[i]);
        if (base == BYTE_MASK) {
            free(encoded_UMI);
            return NULL; // unknown base
        }
        *p |= (base & BASE_MASK) << bit_index;
        bit_index -= BASE_BITS;
        if (bit_index < 0) {
            bit_index = BYTE_SIZE - BASE_BITS;
            p++;
        }
    }
    return encoded_UMI;
}

int main()
{
    char UMI[11] = "ATCGATCGAT";
    printf("length of UMI: %ld\n", strlen(UMI));
    char* encoded_UMI = encode_UMI(UMI);
    printf("%s\n", encoded_UMI);
    return 0;
}