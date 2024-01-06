#include <gtest/gtest.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <SQLiteCpp/SQLiteCpp.h>

using namespace std::this_thread;  // sleep_for, sleep_until
using namespace std::chrono;       // nanoseconds, system_clock, seconds

TEST(DATATWIN, DUCKDB_SQLITE)
{
    // Set the path of the database
    std::string data_dir = "/home/ltruciosr-dev/Documents/astay/datatwin/mining_sim/data/";
    std::string db_file = data_dir + "DbDataTwin.sqlite";
    
    // Open the database in read-only mode.
    SQLite::Database db(db_file, SQLite::OPEN_READONLY);
    
    // SQL query to retrieve mining cuts data.
    std::string query_cuts = "SELECT x,y,z,vertice,period FROM tb_PlanDiario ORDER BY period ASC, vertice ASC";
    SQLite::Statement query(db, query_cuts);

    // SQL query to retrieve mining cuts data.
    while (query.executeStep())
    {
        auto x = query.getColumn(0).getDouble();
        auto y = query.getColumn(1).getDouble();
        auto z = query.getColumn(2).getDouble();
        auto vertice = query.getColumn(3).getInt();
        auto period = query.getColumn(4).getString();
        std::cout << x << "\t" << y << "\t" << z << "\t" << vertice << "\t" << period << std::endl;
    }
    
    // Assert
    ASSERT_EQ(0, 0);
}
