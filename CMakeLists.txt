cmake_minimum_required(VERSION 3.16.0)

project(fastF VERSION 1.0)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Release)

include_directories(/usr/local/lib/R/site-library/Rcpp/include)
include_directories(/usr/local/lib/R/include)
find_package(OpenMP REQUIRED)
# find_package(htslib REQUIRED)
include_directories(./htslib/include)
link_directories(./htslib/lib)

# add_executable(test ./src/v1.1/fastq_filter.c)
# add_executable(fastF ./src/v1.1/main.c ./src/v1.1/fastq_filter.c ./src/v1.1/argparse.c)
# add_executable(fastF ./src/v1.2./main.c ./src/v1.2/fastq_filter.c ./src/v1.2/argparse.c)
add_executable(fastF ./src/v1.2.09/main.c ./src/v1.2.09/fastq_filter.c ./src/v1.2.09/argparse.c ./src/v1.2.09/count.c ./src/v1.2.09/extract.c )
add_executable(test ./test/test.c ./src/v1.2.09/fastq_filter.c ./src/v1.2.09/count.c )
# target_link_libraries(test z pthread)
target_link_libraries(fastF z pthread hts)
target_link_libraries(test z pthread hts)
# target_link_libraries(fastF_openmp z pthread)

# target_link_libraries(fastF PRIVATE OpenMP::OpenMP_C z)
# target_link_libraries(test_omp PRIVATE OpenMP::OpenMP_C)