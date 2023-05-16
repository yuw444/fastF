#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

void *toy(int *b){
    uint8_t *p = (uint8_t *)b;

    for (int i = 0; i < 3; i++) {
        *p = 1 << i;
        p++;
    }
}

int main() {
    printf("Hello, World!\n");

    int *b = (int *)calloc(1, sizeof(int));
    toy(b);

    b = realloc(b, 2*sizeof(int));
    toy(b+1);

    for(int i = 0; i < 2; i++) {
        printf("%d\n", *(b+i));
    }

    free(b);

    printf("Size of bool: %zu bytes\n", sizeof(bool));
    return 0;

    return 0;
}