#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define BASE_BITS 2 // the number of bits to encode a base
#define BYTE_SIZE 8 // the number of bits in a byte
#define BYTE_MASK 0xFF // the mask of a byte
#define BASE_MASK ((1 << BASE_BITS) - 1) // the mask of a base
#define UMI_BYTES_INT (name, umi_length) \
    typedef uint8_t name[(int)((umi_length + 3)/4)]; 
    
uint8_t encode_base(char base) {
    switch (base) {
        case 'A': return 0b00;
        case 'C': return 0b01;
        case 'G': return 0b10;
        case 'T': return 0b11;
        default: return BYTE_MASK; // unknown base
    }
}

uint8_t* compress_sequence(const char* sequence, size_t length, size_t* compressed_length) {
    // calculate the length of compressed sequence
    *compressed_length = (length * BASE_BITS + BYTE_SIZE - 1) / BYTE_SIZE;
    uint8_t* compressed = (uint8_t*)calloc(*compressed_length, sizeof(uint8_t));
    if (compressed == NULL) {
        return NULL;
    }
    uint8_t* p = compressed;
    int bit_index = BYTE_SIZE - BASE_BITS;
    for (size_t i = 0; i < length; i++) {
        uint8_t base = encode_base(sequence[i]);
        if (base == BYTE_MASK) {
            free(compressed);
            return NULL; // unknown base
        }
        *p |= (base & BASE_MASK) << bit_index;
        bit_index -= BASE_BITS;
        if (bit_index < 0) {
            bit_index = BYTE_SIZE - BASE_BITS;
            p++;
        }
    }
    return compressed;
}

typedef uint8_t umi_bytes[3];

bool compress_DNA_seqence(const char* sequence, umi_bytes *compressed) {
    
    if (compressed == NULL) {
        return false;
    }
    int bit_index = BYTE_SIZE - BASE_BITS;
    for(int i = 0; i < 3; i++) {
        compressed[0][i] = 0;
    }
    uint8_t *p = (uint8_t *)compressed;
    for (size_t i = 0; i < strlen(sequence); i++) {
        uint8_t base = encode_base(sequence[i]);
        if (base == BYTE_MASK) {
            return false; // unknown base
        }
        *p |= (base & BASE_MASK) << bit_index;
        bit_index -= BASE_BITS;
        if (bit_index < 0) {
            bit_index = BYTE_SIZE - BASE_BITS;
            p++;
        }
    }
    return true;
}

int main() {
    const char* sequence = "ACGTACGTTT";
    size_t length = strlen(sequence);
    size_t compressed_length;
    uint8_t* compressed = compress_sequence(sequence, length, &compressed_length);
    if (!compressed) {
        printf("Invalid sequence\n");
        return 1;
    }
    printf("Original sequence: %s\n", sequence);
    printf("Compressed sequence: ");
    for (size_t i = 0; i < compressed_length; i++) {
        printf("%02X", compressed[i]);
    }
    printf("\n");
    free(compressed);

    umi_bytes *encoded_umi = calloc(1, sizeof(umi_bytes));

    compress_DNA_seqence("ACGTACGTTT", encoded_umi);

    encoded_umi = realloc(encoded_umi, 2*sizeof(umi_bytes));

    compress_DNA_seqence("AACCGGTTAC", encoded_umi+1);

    printf("encoded_umi: %02X%02X%02X\n", 
        encoded_umi[1][0], 
        encoded_umi[1][1], 
        encoded_umi[1][2]);

    free(encoded_umi);

    return 0;
}
