include(CTest) # Include testing cmake package.
set(TEST_BINARY ${PROJECT_NAME}_test)


add_executable(${TEST_BINARY} 
    test_init.cpp
    test_release1.cpp
    # test_sqlite_duckdb.cpp
    # test_sqlite_sqlitecpp.cpp
    # test_parquet_arrow.cpp
    # test_parquet_duckdb.cpp
)

target_link_libraries(${TEST_BINARY}
    gtest
    gtest_main
    SQLiteCpp
    sqlite3
    duckdb
    pthread
    mining_sim
    ) 

add_test(
    NAME ${TEST_BINARY}
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/${TEST_BINARY})