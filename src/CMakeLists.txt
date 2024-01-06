add_library(mining_sim STATIC mining_sim.cpp)

target_link_libraries (mining_sim 
    SQLiteCpp
    sqlite3
    pthread
    duckdb
    )