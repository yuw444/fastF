#include "utils.h"

size_t *GetSeqInt(size_t start, size_t end, size_t step)
{
  if (step == 0)
  {
    printf("Step cannot be 0");
    exit(1);
  }

  if ((end - start) * step < 0)
  {
    printf("Step must be positive when start < end, or negative when start > end");
    exit(1);
  }

  size_t *seq = (size_t *)calloc((size_t)(end + step -1 - start)/step + 1, sizeof(size_t));

  size_t i = 0;
  for (size_t j = start; j <= end; j += step)
  {
    seq[i] = j;
    i++;
  }

  return seq;
}

size_t *SampleInt(size_t *arrayIn, size_t nTotal, size_t nSample, unsigned int replace, unsigned int seed)
{

  init_genrand(seed);

  size_t *sampleOut = calloc(nSample, sizeof(size_t));

  if (replace == 0)
  {
    if (nSample > nTotal)
    {
      free(sampleOut);
      printf("Sample size must be smaller than population size when sampling without replacement.");
      exit(1);
    }

    size_t *arrayInCopy = calloc(nTotal, sizeof(size_t));
    memcpy(arrayInCopy, arrayIn, nTotal * sizeof(size_t));

    if (nTotal == nSample) {
        free(sampleOut);
        return arrayInCopy;
    }

    for (size_t i = 0; i < nSample; i++)
    {
      size_t index = genrand_int32() % nTotal;
      sampleOut[i] = arrayInCopy[index];
      if (index != (nTotal - 1))
      {
        arrayInCopy[index] = arrayInCopy[nTotal - 1];
      }
      nTotal--;
    }

    free(arrayInCopy);
  }
  else
  {
    for (size_t i = 0; i < nSample; i++)
    {
      sampleOut[i] = arrayIn[genrand_int32() % nTotal];
    }
  }

  return sampleOut;
}

void PrintArrayInt(size_t *array, size_t n)
{
  for (size_t i = 0; i < n; i++)
  {
    printf("%d ", array[i]);
  }
  printf("\n");
}

int vsI(const void *a, const void *b)
{
  return (*(size_t *)a - *(size_t *)b);
}