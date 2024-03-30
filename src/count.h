#ifndef COUNT_H
#define COUNT_H
#include "extract.h"

// get both cell barcode from R1
node *cell_counts (gzFile R1_file, size_t len_elements, size_t len_umi);

#endif