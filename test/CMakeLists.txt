set(TEST_SOURCES_UTILS
    test_data.c
    ../src/bamonly/data.h
)

set(CONVERT
    # ./convert_bam.c
    ./hash_cellbarcode.c
    ../external/hashtable/hashtable.c
)

set(FASTF
    
    ./fastF.c
    ../external/hashtable/hashtable.c)


set(CODING
    
    ./test_decode.c
    ../src/v1.2.0999/bam2db_pure.c
    ../external/hashtable/hashtable.c
    )

set(TABLE2GZ

    ./test_sql_table2gz.c)

set(CDTABLE

    ./test_sql_create_drop_table.c)

add_executable(test_data ${TEST_SOURCES_UTILS})
add_executable(convert ${CONVERT})
add_executable(fastf ${FASTF})
add_executable(coding ${CODING})
add_executable(table2gz ${TABLE2GZ})
add_executable(cdtable ${CDTABLE})
target_link_libraries(test_data criterion m)
target_link_libraries(convert hts z criterion)
target_link_libraries(fastf hts z sqlite3 criterion)
target_link_libraries(coding hts z sqlite3 criterion)
target_link_libraries(table2gz z sqlite3)
target_link_libraries(cdtable z sqlite3)

# add_test(NAME MyTests COMMAND test_data)

add_custom_target(TestData
    COMMAND valgrind  --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes  --log-file=vg_logfile_%p.out ./test_data
    DEPENDS test_data
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(Convert
    COMMAND ./convert
    DEPENDS convert
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(FastF
    COMMAND ./fastf
    DEPENDS fastf
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(Coding
    # COMMAND valgrind --show-leak-kinds=all --trace-children=yes --leak-check=full -v --track-origins=yes --dsymutil=yes  --log-file=vg_logfile_%p.out ./coding
    COMMAND ./coding
    DEPENDS coding
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
