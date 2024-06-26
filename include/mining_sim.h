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
#include "geometry/block.h"
#include "geometry/drawer.h"
#include "geometry/polyline.h"
#include "planner/plan.h"
#include "mining_sim.h"

namespace bg = boost::geometry;

namespace astay
{
    class Datatwin
    {
    private:
        typedef double float_t;
        typedef bg::model::d2::point_xy<float_t> point_t;
        typedef bg::model::polygon<point_t> polygon_t;
        typedef bg::model::box<point_t> box_t;
        typedef bg::model::linestring<point_t> linestring_t;
        std::string db_filepath_;
        std::string mb_filepath_;
        std::string tp_filepath_;
        std::string svg_dir_;
        float_t xOffset_, yOffset_;
        int mb_height_{10};
        std::vector<MiningCut> cuts_;
        std::vector<MiningBlock> blocks_;
        std::vector<Polyline> polylines_;
        std::vector<MiningGeoPoly> geopolygons_;
        std::unordered_map<int, std::vector<int>> cut_id_by_level_;
        std::unordered_map<int, std::vector<int>> polyline_id_by_level_;
        std::unordered_map<int, std::vector<int>> geopoly_id_by_cut_id_;
        std::unordered_map<int, std::vector<int>> block_id_by_cut_id_;
        // std::unordered_map<int, std::vector<int>> geopoly_id_by_level_;
        // std::map<int, std::vector<MiningBlock>> blocks_by_level_;
        // std::map<std::string, std::vector<MiningCut>> cuts_by_loader_;

    public:
        Datatwin();
        /* Set DB parameters */
        void setDBfilepath(const std::string &filepath);
        void setModelBlockfilepath(const std::string &filepath);
        void setTopographyfilepath(const std::string &filepath);
        void setSVGDir(const std::string &filepath);
        /* Read data */
        void ReadMiningCuts(bool is_verbose = false);
        void ReadMiningGeoPolygons(bool is_verbose = false);
        void ReadTopography(bool is_verbose = false);
        void ReadMiningBlocks(bool is_verbose = false);
        // void ReadPlan(bool is_verbose = false);
        /* Compute operations */
        void IntersectMiningCutsWithPolylines(bool is_verbose = false);
        void Draw(bool is_verbose = false);

    private:
        template <class PolygonType>
        PolygonType ApplyOffsetToPolygon(const PolygonType &polygon, point_t &centroid)
        {
            PolygonType out_polygon;
            polygon_t poly = polygon.geometry();
            for (auto it = boost::begin(bg::exterior_ring(poly));
                 it != boost::end(bg::exterior_ring(poly));
                 ++it)
            {
                float_t x = bg::get<0>(*it) - centroid.x();
                float_t y = bg::get<1>(*it) - centroid.y();
                point_t p{x, y};
                out_polygon.add(p);
            }
            out_polygon.correct();
            return out_polygon;
        }
        Polyline ApplyOffsetToPolyline(const Polyline &polyline, point_t &centroid)
        {
            Polyline out_polyline;
            linestring_t linestring = polyline.linestring();
            for (auto it = boost::begin(linestring); it != boost::end(linestring); ++it)
            {
                float_t x = bg::get<0>(*it) - centroid.x();
                float_t y = bg::get<1>(*it) - centroid.y();
                point_t p{x, y};
                out_polyline.add(p);
            }
            return out_polyline;
        }
        void finalizeAndStoreCut(MiningCut& cut, int level, int day);
        void finalizeAndStorePolyline(Polyline& polyline, int level);
        void finalizeAndStoreGeoPoly(MiningGeoPoly& geopoly, int level, std::string& name);
        bool isLeftOfLine(point_t& line_a, point_t& line_b, point_t& point_c);
        float_t distanceFromLine(point_t &line_a, point_t &line_b, point_t &point_c);
        box_t computeBoundingArea(float margin = 1.5, bool is_verbose = false);
        linestring_t linestringFromPolyline(Polyline& polyline, point_t point_a, point_t point_b);
    };
}