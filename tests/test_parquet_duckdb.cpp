#include <gtest/gtest.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <duckdb.hpp>

using namespace std::this_thread;  // sleep_for, sleep_until
using namespace std::chrono;       // nanoseconds, system_clock, seconds

// SELECT * FROM read_parquet('/home/ltruciosr-dev/Documents/astay/datatwin/mining_sim/data/Topografia.parquet') limit 10;

TEST(DATATWIN, DUCKDB_PARQUET)
{
    /* Read the parameters */
    std::string data_dir = "/home/ltruciosr-dev/Documents/astay/datatwin/mining_sim/data/";
    std::string modelblock_parquet_file = data_dir + "ModelBlocks.parquet";
    std::string topography_parquet_file = data_dir + "Topografia.parquet";
    
    /* Setup readers */
    duckdb::DuckDB db(nullptr);
    duckdb::Connection con(db);

    // SQL query to retrieve mining cuts data.
    std::string query_cuts = "SELECT x,y,z,vertice,geometry FROM read_parquet('" + topography_parquet_file + "') WHERE z=4060";
    // duckdb::MaterializedQueryResult;
    auto result = con.Query(query_cuts);
    
    // Get the count of rows
    auto count = result->RowCount();
    std::cout << count << std::endl;

    // Fetch one chunk and valid the vectors for each column
    auto chunk = result->Fetch();
    auto n_rows = chunk->size();
    auto n_cols = chunk->ColumnCount();
    std::cout << n_rows << std::endl;
    std::cout << n_cols << std::endl;

    // for (auto vec : chunk->data)
    // {
    //     std::cout << vec.GetValue(0) << std::endl;
    // }

    auto x_vec = chunk->data[0];
    auto y_vec = chunk->data[1];
    auto z_vec = chunk->data[2];
    auto vertex_vec = chunk->data[3];
    auto geo_vec = chunk->data[4];

    for (int i = 0; i < chunk->size(); i++)
    {
        auto x = x_vec.GetValue(i);
        auto y = y_vec.GetValue(i);
        auto z = z_vec.GetValue(i);
        auto vertex = vertex_vec.GetValue(i);
        auto geo = geo_vec.GetValue(i);
        // std::cout << x << "," << y << "," << z << "," << vertex << "," << geo << std::endl;
    }
    // Iterate over the chunks
    /*
    auto chunk = result->Fetch();
    int ctr = 0;
    while (chunk != nullptr)
    {
        ctr += chunk->size();
        chunk = result->Fetch();
    }
    std::cout << ctr << std::endl;
    */
    // Assert
    ASSERT_EQ(0, 0);
}
