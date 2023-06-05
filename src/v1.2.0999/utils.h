#ifndef UTILS_H
#define UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mt19937ar.h"

size_t *GetSeqInt(size_t start, size_t end, size_t step);

size_t *SampleInt(size_t *arrayIn, size_t nTotal, size_t nSample, unsigned int replace, unsigned int seed);

void PrintArrayInt(size_t *array, size_t n);

int vsI(const void *a, const void *b);
#endif