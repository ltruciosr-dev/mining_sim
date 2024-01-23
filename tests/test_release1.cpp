#include <gtest/gtest.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <filesystem>

#include "mining_sim.h"

using namespace std::this_thread;  // sleep_for, sleep_until
using namespace std::chrono;       // nanoseconds, system_clock, seconds

TEST(DATATWIN, RELEASE_1)
{
    /* Read the parameters */
    std::string data_dir = "/home/ltruciosr-dev/Documents/astay/datatwin/mining_sim/data/";
    std::string db_file = data_dir + "DbDataTwin.sqlite";
    std::string modelblock_parquet_file = data_dir + "ModelBlocks.parquet";
    std::string topography_parquet_file = data_dir + "Topografia.parquet";
    std::string svg_dir = data_dir + "svg/release_1/";
    std::cout << "RELEASE 1 | Store SVG on -> " << svg_dir << std::endl;
    
    /* Setup parameters */
    astay::Datatwin datatwin;
    datatwin.setDBfilepath(db_file);
    datatwin.setTopographyfilepath(topography_parquet_file);
    datatwin.setModelBlockfilepath(modelblock_parquet_file);
    datatwin.setSVGDir(svg_dir);
    // datatwin.setModelBlocksDB(modelblocks_db);
    
    /* Read Mining Polygons and Topography */
    datatwin.ReadMiningCuts(false);
    datatwin.ReadTopography(false);
    datatwin.ReadMiningGeoPolygons(false);
    // datatwin.ReadMiningBlocks(false);
    
    /* Process Information */
    datatwin.IntersectMiningCutsWithPolylines(true);
    // -- compute
    // datatwin.ComputeSlices(slice_dist, kShowVerbose);
    // datatwin.ComputeFingerIntersections(kDontVerbose);
    // -- draw
    datatwin.Draw(false);

    // Assert
    ASSERT_EQ(0, 0);
}
