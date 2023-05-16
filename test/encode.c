#include "./encode.h"

// ENCODED_BYTE_INT(umi_bytes, 10);
// bool compress_UMI(const char *sequence, umi_bytes *compressed) 
// { 
//     if (compressed == ((void *)0)){ return 0; } 
//     int bit_index = 8 - 2; 
//     for (int i = 0; i < 3; i++) 
//     { 
//         compressed[0][i] = 0; 
//     } 
//     uint8_t *p = (uint8_t *)compressed; 
//     for (size_t i = 0; i < strlen(sequence); i++) 
//     { 
//         uint8_t base = (sequence[i] == 'A' ? 0b00 : sequence[i] == 'C' ? 0b01 : sequence[i] == 'G' ? 0b10 : sequence[i] == 'T' ? 0b11 : 0xFF); 
//         if (base == 0xFF) { return 0; } 
//         *p |= (base & ((1 << 2) - 1)) << bit_index; 
//         bit_index -= 2; 
//         if (bit_index < 0) 
//         { 
//             bit_index = 8 - 2; 
//             p++; 
//         } 
//     } 

//     return 1; 
// }

int main() {
    const char* sequence = "ACGTACGTTT";
    ENCODED_BYTE_INT(umi_bytes, 10);
    ENCODED_FUNC_INT(umi_bytes, compress_UMI);

    umi_bytes *compressed = calloc(1, sizeof(umi_bytes));

    if(compress_UMI("AACCTTGGTT", compressed))
        printf("%02X %02X %02X\n", compressed[0][0], compressed[0][1], compressed[0][2]);


    compressed = realloc(compressed, 2*sizeof(umi_bytes));

    compress_UMI("AACTTTGGTT", compressed+1);

    printf("%02X %02X %02X\n", compressed[1][0], compressed[1][1], compressed[1][2]);    

    printf("size of (size_t) = %lu\n", sizeof(size_t));
    printf("size of (unsigned int) = %lu\n", sizeof(unsigned int));

    free(compressed);

    return 0;
}
