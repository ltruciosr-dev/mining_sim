cmake_minimum_required(VERSION 3.1)
project(mining_sim VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

find_package(Threads)
find_package(DuckDB)
find_package(SQLite3)

include_directories(
    ${PROJECT_SOURCE_DIR}/thirdparty/SQLiteCpp/include
    ${PROJECT_SOURCE_DIR}/include
    ${SQLite3_INCLUDE_DIRS}
    ${DuckDB_INCLUDE_DIRS}
)

set(SQLITECPP_RUN_CPPCHECK OFF CACHE BOOL "" FORCE)
set(SQLITECPP_RUN_CPPLINT OFF CACHE BOOL "" FORCE)
set(SQLITECPP_USE_STATIC_RUNTIME OFF CACHE BOOL "" FORCE)

add_subdirectory(thirdparty/SQLiteCpp)
add_subdirectory(src)
add_subdirectory(tests)
