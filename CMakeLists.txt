cmake_minimum_required(VERSION 3.16.0)

project(fastF VERSION 1.0)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

include_directories(/usr/local/lib/R/site-library/Rcpp/include)
include_directories(/usr/local/lib/R/include)
# find_package(OpenMP REQUIRED)

add_executable(test ./src/v1.0.9/fastq_filter.c)
add_executable(test_parse ./src/v1.0.9/fastq_parse.c ./src/v1.0.9/argparse.c)

target_link_libraries(test z pthread)

# target_link_libraries(test_forest PRIVATE OpenMP::OpenMP_C m)
# target_link_libraries(test_omp PRIVATE OpenMP::OpenMP_C)