#include <gtest/gtest.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <duckdb.hpp>

using namespace std::this_thread;  // sleep_for, sleep_until
using namespace std::chrono;       // nanoseconds, system_clock, seconds

TEST(DATATWIN, DUCKDB_SQLITE)
{
    /* Read the parameters */
    std::string data_dir = "/home/ltruciosr-dev/Documents/astay/datatwin/mining_sim/data/";
    std::string db_file = data_dir + "DbDataTwin.sqlite";
    
    /* Setup readers */
    duckdb::DuckDB db(db_file);
    duckdb::Connection con(db);

    // SQL query to retrieve mining cuts data.
    std::string query_cuts = "SELECT x,y,z,vertice,period FROM tb_PlanDiario ORDER BY period ASC, vertice ASC";
    auto result = con.Query(query_cuts);
    result->Print();
    
    // Assert
    ASSERT_EQ(0, 0);
}
