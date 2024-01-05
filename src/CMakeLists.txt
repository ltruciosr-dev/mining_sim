add_library(kaz_tools STATIC kaz_mining.cpp)
add_library(kaz_sort STATIC kaz_sort.cpp)

target_link_libraries (kaz_tools ${SQLite3_LIBRARIES})
target_link_libraries (kaz_sort ${SQLite3_LIBRARIES})