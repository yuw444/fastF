#include <stdio.h>
#include <stdlib.h>

#define OPTION 0x16 // 00010000

int main()
{

    for (int i = 0; i < 64; i++)
    {
        if (!(i & (1 << 4)))
            printf("%d, SET\n", i);
        else
            printf("%d, NOT SET\n", i);
    }

    return 0;
}