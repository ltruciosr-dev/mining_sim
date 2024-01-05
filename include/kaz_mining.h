#pragma once

#include <string>
#include <iostream>
#include <sqlite3.h>
#include <map>

#include "geometry/polygon.h"
#include "geometry/drawer.h"
#include "geometry/finger.h"
#include "geometry/model_block.h"
#include "SQLiteCpp/SQLiteCpp.h"

#include "kaz_mining.h"

namespace astay
{
    class KazMining
    {
    private:
        typedef double float_t;
        typedef bg::model::d2::point_xy<float_t> point_xy;
        typedef bg::model::polygon<point_xy> polygon_t;
        SQLite::Database db_;
        SQLite::Transaction transaction_;
        // -- parameters
        std::vector<std::string> loaders_{"EX201", "EX202", "EX203", "EX204"};
        int project_id_;
        std::string svg_dir_;
        std::string project_db_;
        std::string modelblocks_db_;
        float_t x_offset_, y_offset_;
        std::map<std::string, std::vector<Finger>> fingers_by_loader_;
        std::map<std::string, std::vector<int>> fingers_order_by_loader_;
        // -- compute
        void ComputeOreIntersections(const polygon_t &cut_inter, int order, SQLite::Database &db);
        void ComputeBlockIntersections(const polygon_t &ore_inter, int order, SQLite::Database &db);
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
        // -- output
        float_t t_area_;
        float_t t_cu_;
        float_t t_ton_;
        int num_slice_;
        int finger_id_;
        std::string loader_;
        std::string ore_name_;
        std::array<float_t, 3> ore_position_;

    public:
        KazMining(const std::string &db_name);
        // -- set
        void setProjectID(int project_id);
        void setProjectDB(const std::string &db_name);
        void setModelBlocksDB(const std::string &db_name);
        void setSVGDir(const std::string &svg_dir);
        // -- read
        void ReadCuts(bool verbose = false);
        void ReadOres(bool verbose = false);
        void ReadBlocks(bool verbose = false);
        // -- compute
        void ComputeSlices(float_t slice_dist, bool verbose = false);
        void ComputeFingerIntersections(bool verbose = false);
        void Draw(bool verbose = false);
        };
}