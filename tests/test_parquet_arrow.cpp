#include <gtest/gtest.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <arrow/io/api.h>
#include <parquet/stream_reader.h>

using namespace std::this_thread;  // sleep_for, sleep_until
using namespace std::chrono;       // nanoseconds, system_clock, seconds

TEST(DATATWIN, ARROW_PARQUET)
{
    /* Read the parameters */
    std::string data_dir = "/home/ltruciosr-dev/Documents/astay/datatwin/mining_sim/data/";
    std::string modelblock_parquet_file = data_dir + "ModelBlocks.parquet";
    std::string topography_parquet_file = data_dir + "Topografia.parquet";
    
    /* Setup parameters */
    std::shared_ptr<arrow::io::ReadableFile> infile;
    PARQUET_ASSIGN_OR_THROW(
      infile,
      arrow::io::ReadableFile::Open(topography_parquet_file));

    parquet::StreamReader stream{parquet::ParquetFileReader::Open(infile)};

    double x;
    double y;
    double z;
    int64_t vertice;
    std::string type_obj;

    std::vector<double> x_curves, y_curves, z_curves;
    std::vector<int64_t> vertex_numbers;

    while ( !stream.eof() )
    {
       stream >> x >> y >> z >> vertice >> type_obj >> parquet::EndRow;
       if (z = 4059)
       {
            x_curves.push_back(x);
            y_curves.push_back(y);
            z_curves.push_back(z);
            vertex_numbers.push_back(vertice);
       }
    }
    std::cout << x << "," << y << "," << z << "," << vertice << "," << type_obj << std::endl;
    // Assert
    ASSERT_EQ(0, 0);
}
