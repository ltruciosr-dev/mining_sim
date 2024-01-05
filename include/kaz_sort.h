#pragma once

#include <string>
#include <iostream>
#include <sqlite3.h>
#include <map>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include "SQLiteCpp/Database.h"
#include "SQLiteCpp/Transaction.h"

#include "geometry/polygon.h"
#include "geometry/model_block.h"
#include "geometry/drawer.h"
#include "planner/plan.h"
#include "kaz_sort.h"

namespace bg = boost::geometry;

namespace astay
{
    class KazSort
    {
    private:
        typedef double float_t;
        typedef bg::model::d2::point_xy<float_t> point_xy;
        typedef bg::model::polygon<point_xy> polygon_t;
        SQLite::Database db_;
        SQLite::Transaction transaction_;
        // -- parameters
        int project_id_;
        bool is_sorted_;
        std::string project_db_;
        std::string modelblocks_db_;
        std::string svg_dir_;
        float_t x_offset_, y_offset_;
        MiningPlan plan_;
        std::vector<int> row_ids_;
        std::vector<int> row_cuts_;
        std::vector<MiningCut> cuts_;
        std::vector<std::string> row_loaders_;
        std::map<int, std::vector<MiningOre>> ores_by_level_;
        std::map<int, std::vector<MiningBlock>> blocks_by_level_;
        std::map<std::string, std::vector<MiningCut>> cuts_by_loader_;
        // -- compute
        template <class PolygonType>
        PolygonType FixCoordinates(const PolygonType &polygon, point_xy &centroid)
        {
            PolygonType out_polygon;
            polygon_t poly = polygon.geometry();
            for (auto it = boost::begin(bg::exterior_ring(poly));
                 it != boost::end(bg::exterior_ring(poly));
                 ++it)
            {
                float_t x = bg::get<0>(*it) - centroid.x();
                float_t y = bg::get<1>(*it) - centroid.y();
                point_xy p{x, y};
                out_polygon.add(p);
            }
            out_polygon.correct();
            return out_polygon;
        }
        MiningBlock FixCoordinates(const MiningBlock &block, point_xy &centroid)
        {
            MiningBlock block_moved;
            float_t block_xmin = bg::get<bg::min_corner, 0>(block.geometry());
            float_t block_ymin = bg::get<bg::min_corner, 1>(block.geometry());
            float_t block_xmax = bg::get<bg::max_corner, 0>(block.geometry());
            float_t block_ymax = bg::get<bg::max_corner, 1>(block.geometry());
            point_xy p_min{block_xmin - centroid.x(), block_ymin - centroid.y()};
            point_xy p_max{block_xmax - centroid.x(), block_ymax - centroid.y()};
            block_moved.corners(p_min, p_max);
            block_moved.set_id(block.id());
            block_moved.correct();
            return block_moved;
        }

    public:
        KazSort(std::string db_name);
        // -- set
        void setProjectID(int project_id);
        void setProjectDB(const std::string &db_name);
        void setModelBlocksDB(const std::string &db_name);
        void setSVGDir(const std::string &svg_dir);
        // -- read
        void ReadCuts(bool verbose = false);
        void ReadOres(bool verbose = false);
        void ReadBlocks(bool verbose = false);
        void ReadPlan(bool verbose = false);
        // -- compute
        void SortCuts(bool verbose = false);
        void ValidateDB(bool verbose = false);
        void PopulateDB(bool verbose = false);
        void Draw(bool verbose = false);
    };
}