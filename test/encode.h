
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define COLOR_RED "\x1b[31m"
#define COLOR_GREEN "\x1b[32m"
#define COLOR_YELLOW "\x1b[33m"
#define COLOR_BLUE "\x1b[34m"
#define COLOR_RESET "\x1b[0m"
#define BASE_BITS 2                      // the number of bits to encode a base
#define BYTE_SIZE 8                      // the number of bits in a byte
#define BYTE_MASK 0xFF                   // the mask of a byte
#define BASE_MASK ((1 << BASE_BITS) - 1) // the mask of a base
#define ENCODE_BASE(base) (base == 'A' ? 0b00 : base == 'C' ? 0b01 : base == 'G' ? 0b10 : base == 'T' ? 0b11 : BYTE_MASK)

// define a new data type to store the encoded sequence with fixed length (length is the number of bases)
#define ENCODED_BYTE_INT(encoded_byte_type_name, seq_length) \
    typedef uint8_t encoded_byte_type_name[(int)((seq_length + BYTE_SIZE/BASE_BITS - 1 ) / (BYTE_SIZE/BASE_BITS))];

// define a function to encode a sequence to any data type
// this function encode a seqence to the data type that sizeof(data_type) greater than or equal to strlen(sequence)/(BYTE_SIZE/BASE_BITS)
// The ideal data type is uint8_t, uint16_t, uint32_t, uint64_t or their array
// for example 
// uint8_t is the ideal data type for a sequence with length <= 4
// uint16_t is the ideal data type for a sequence with length <= 8
// uint32_t is the ideal data type for a sequence with length <= 16
// uint64_t is the ideal data type for a sequence with length <= 32
#define ENCODED_FUNC_INT(encoded_byte_type, function_name) \
bool function_name(const char *sequence, encoded_byte_type *compressed) {\
    if (compressed == NULL){\
        return false; \
    }\
    if (strlen(sequence) > sizeof(encoded_byte_type)*BYTE_SIZE/BASE_BITS) {\
        printf(COLOR_RED "ERROR: The encoded sequence will be longer than size of " COLOR_GREEN #encoded_byte_type COLOR_RED ", check how " COLOR_GREEN #encoded_byte_type COLOR_RED " is defined\n" COLOR_RESET);\
        return false; \
    }\
    int bit_index = BYTE_SIZE - BASE_BITS; \
    memset(compressed, 0, sizeof(encoded_byte_type));\
    uint8_t *p = (uint8_t *)compressed;\
    for (size_t i = 0; i < strlen(sequence); i++) { \
        uint8_t base = ENCODE_BASE(sequence[i]); \
        if (base == BYTE_MASK) { \
            return false; \
        } \
        *p |= (base & BASE_MASK) << bit_index; \
        bit_index -= BASE_BITS; \
        if (bit_index < 0) {\
            bit_index = BYTE_SIZE - BASE_BITS; \
            p++; \
        } \
    } \
    return true; \
}
