#include "../src/utils.h"

int main()
{
    size_t *seq = GetSeqInt(0, 100, 1);
    PrintArrayInt(seq, 101);
    size_t *sample = SampleInt(seq, 101, 101, 0, 1);
    qsort(sample, 20, sizeof(size_t), vsI);
    PrintArrayInt(sample, 101);

    printf("Size of int: %d\n", sizeof(int));

    free(seq);
    free(sample);
    return 0;
}
