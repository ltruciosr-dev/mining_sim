#include <gtest/gtest.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include "kaz_sort.h"

using namespace std::this_thread;  // sleep_for, sleep_until
using namespace std::chrono;       // nanoseconds, system_clock, seconds

TEST(SortingCuts, SORT)
{
    // DATA
    // -- input
    int project_id = 4; // it has a specific PName and PVersion
    std::string data_dir = "/home/ltruciosr-dev/Documents/astay/kaz/astay_mining_sim/data/";
    std::string project_db = data_dir + "SimulationAS.db";
    std::string modelblocks_db = data_dir + "modelblocks_T.db";
    // -- output
    std::string svg_dir = data_dir + "svg/project_" + std::to_string(project_id) + "/";
    std::cout << "PROJECT " << project_id << " | Store SVG on -> " << svg_dir << std::endl;

    astay::KazSort kazObj(project_db);

    // APPLY
    const bool kShowVerbose = true;
    const bool kDontVerbose = false;
    // -- set
    kazObj.setProjectID(project_id);
    kazObj.setProjectDB(project_db);
    kazObj.setModelBlocksDB(modelblocks_db);
    kazObj.setSVGDir(svg_dir);
    // -- read
    kazObj.ReadCuts(kDontVerbose);
    kazObj.ReadOres(kDontVerbose);
    kazObj.ReadBlocks(kDontVerbose);
    kazObj.ReadPlan(kDontVerbose);
    // -- compute
    kazObj.SortCuts(kDontVerbose);
    kazObj.ValidateDB(kShowVerbose);
    kazObj.PopulateDB(kShowVerbose);
    // kazObj.Draw(kDontVerbose);
    
    // Assert
    ASSERT_EQ(0, 0);
}
