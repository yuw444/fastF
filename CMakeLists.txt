cmake_minimum_required(VERSION 3.16.0)

project(fastF VERSION 1.0)

set(CMAKE_C_STANDARD 99)
set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Release)

find_package(OpenMP REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/external/htslib/include)
link_directories(${CMAKE_SOURCE_DIR}/external/htslib/lib)
include_directories(${CMAKE_SOURCE_DIR}/external/criterion/include)
link_directories(${CMAKE_SOURCE_DIR}/external/criterion/lib)
include_directories(${CMAKE_SOURCE_DIR}/external/argparse/)
include_directories(${CMAKE_SOURCE_DIR}/external/mt19937ar.sep)
include_directories(${CMAKE_SOURCE_DIR}/external/hashtable)


set(SOURCES
    ./src/v1.2.0999/main.c
    ./src/v1.2.0999/filter.c
    external/argparse/argparse.c
    external/mt19937ar.sep/mt19937ar.c 
    ./src/v1.2.0999/count.c
    ./src/v1.2.0999/extract.c
    external/hashtable/hashtable.c
    ./src/v1.2.0999/bam2db_ds.c
)

add_executable(fastF ${SOURCES})

target_link_libraries(fastF z pthread hts sqlite3)

add_subdirectory(./test)