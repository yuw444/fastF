#ifndef COUNT_H
#define COUNT_H

#include "fastq_filter.h"

node *cell_counts (gzFile R1_file, size_t len_cellbarcode, size_t len_umi);

#endif