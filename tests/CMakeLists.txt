# add_subdirectory(/usr/lib/
#     ${PROJECT_BINARY_DIR}/gtest)
include(CTest) # Include testing cmake package.
set(TEST_BINARY ${PROJECT_NAME}_test)


add_executable(${TEST_BINARY} 
    test_init.cpp
    # test_sort.cpp
    test_mining.cpp
)

target_link_libraries(${TEST_BINARY}
    gtest gtest_main SQLiteCpp sqlite3 pthread kaz_tools kaz_sort) 

add_test(
    NAME ${TEST_BINARY}
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/${TEST_BINARY})