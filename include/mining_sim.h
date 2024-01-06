/*
 * mining_sim.h
 *
 * URL:      https://github.com/ltruciosr-dev/mining_sim
 * Version:  1.0
 *
 * Copyright (C) 2024 Luis Trucios
 * All rights reserved.
 *
 * mining_sim is distributed under the BSD 3-Clause license, see LICENSE for details.
 *
 */
#pragma once

#include <string>
#include <iostream>
#include <sqlite3.h>
#include <map>
#include <unordered_map>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include "SQLiteCpp/Database.h"
#include "SQLiteCpp/Transaction.h"
/* Libraries developed by our team */
#include "geometry/polygon.h"
#include "geometry/model_block.h"
#include "geometry/drawer.h"
#include "planner/plan.h"
#include "mining_sim.h"

namespace bg = boost::geometry;

namespace astay
{
    class Datatwin
    {
    private:
        typedef double float_t;
        typedef bg::model::d2::point_xy<float_t> point_xy;
        typedef bg::model::polygon<point_xy> polygon_t;
        typedef bg::model::box<point_xy> box_t;
        std::string db_filepath_;
        std::string mb_filepath_;
        std::string tp_filepath_;
        // std::string svg_dir_;
        float_t xOffset_, yOffset_;
        float_t mb_height_;
        // MiningPlan plan_;
        // std::vector<int> row_ids_;
        // std::vector<int> row_cuts_;
        std::vector<MiningCut> cuts_;
        std::vector<MiningGeoPoly> geopolygons_;
        // std::vector<std::string> row_loaders_;
        std::unordered_map<int, std::vector<int>> cut_id_by_level_;
        std::unordered_map<int, std::vector<int>> geopoly_id_by_cut_id_;
        // std::unordered_map<int, std::vector<int>> geopoly_id_by_level_;
        // std::map<int, std::vector<MiningBlock>> blocks_by_level_;
        // std::map<std::string, std::vector<MiningCut>> cuts_by_loader_;
        /* fix coordinates for bf::model::polygon inheritance */

    public:
        Datatwin();
        /* Set DB parameters */
        // void setProjectID(int project_id);
        void setDBfilepath(const std::string &filepath);
        void setModelBlockfilepath(const std::string &filepath);
        void setTopographyfilepath(const std::string &filepath);
        // void setModelBlocksDB(const std::string &db_name);
        // void setSVGDir(const std::string &svg_dir);
        /* Read data */
        void ReadMiningCuts(bool is_verbose = false);
        void ReadMiningGeoPolygons(bool is_verbose = false);
        void ReadMiningBlocks(bool is_verbose = false);
        // void ReadTopography(bool is_verbose = false);
        // void ReadPlan(bool is_verbose = false);
        /* Compute operations */
        // void SortCuts(bool is_verbose = false);
        // void ValidateDB(bool is_verbose = false);
        // void PopulateDB(bool is_verbose = false);
        // void Draw(bool is_verbose = false);
    
    private:
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
            block_moved.setId(block.id());
            block_moved.correct();
            return block_moved;
        }
        void finalizeAndStoreCut(MiningCut& cut, int level, int day);
        void finalizeAndStoreGeoPoly(MiningGeoPoly& geopoly, int level, std::string& name);
        box_t computeBoundingArea(float margin = 1.5, bool is_verbose = false);
    };
}