#pragma once
#include <iostream>
#include <fstream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include "geometry/block.h"
#include "geometry/polygon.h"
#include "geometry/polyline.h"
// #include "geometry/finger.h"

namespace bg = boost::geometry;

namespace astay
{
    class Drawer
    {
    private:
        typedef double float_t;
        typedef bg::model::d2::point_xy<float_t> point_t;
        typedef bg::model::polygon<point_t> polygon_t;
        typedef bg::model::box<point_t> box_t;
        std::ofstream svg_;
        bg::svg_mapper<point_t> mapper_;
        int width_, height_;

    public:
        Drawer(const std::string &svg_file, int width = 1600, int height = 1600, std::string const &width_height = "width=\"100%\" height=\"100%\"")
            : svg_(svg_file),
              mapper_(svg_, width, height, width_height)
        {
            width_ = width;
            height_ = height;
            AddFrame();
        };
        void AddFrame()
        {
            point_t min_corner(-width_ / 2, -height_ / 2);
            point_t max_corner(width_ / 2, height_ / 2);
            box_t outer_box(min_corner, max_corner);
            mapper_.add(outer_box);
            mapper_.map(outer_box, "opacity:0.5;fill:none;stroke:rgb(0,0,0);stroke-width:2");
        }
        void AddPoint(point_t point, int radius = 1)
        {
            mapper_.add(point);
            mapper_.map(point, "fill-opacity:1;fill:#02557A;stroke:#012231;stroke-width:0.5", 5);
        }
        void AddNumber(int number, point_t position)
        {
            std::string text = std::to_string(number);
            double offset_x{0}, offset_y{8};
            if (number < 10)
                offset_x = -6;
            else
                offset_x = -12;
            mapper_.text(position, text, "fill-opacity:1;fill:#002255;stroke:#000000;stroke-width:0.2;font-size:10px", offset_x, offset_y);
        }
        void AddCut(const MiningCut &cut, int extraction_day = -1)
        {
            mapper_.add(cut.geometry());
            switch (extraction_day)
            {
            case 1: // verde
                mapper_.map(cut.geometry(), "fill-opacity:0.1;fill:#00FF00;stroke:#00FF00;stroke-width:0.25");
                break;
            case 2: // azul claro
                mapper_.map(cut.geometry(), "fill-opacity:0.1;fill:#0000FF;stroke:#0000FF;stroke-width:0.25");
                break;
            case 3: // anaranjado
                mapper_.map(cut.geometry(), "fill-opacity:0.1;fill:#FF0000;stroke:#FF0000;stroke-width:0.25");
                break;
            default:
                mapper_.map(cut.geometry(), "fill-opacity:0.1;fill:#00008B;stroke:#00008B;stroke-width:0.25");
                break;
            }
        }
        void AddPolyline(const Polyline &polyline)
        {
            mapper_.add(polyline.linestring());
            mapper_.map(polyline.linestring(), "opacity:0.85;fill:#51D1F6;stroke:#000000;stroke-width:0.10");
        }
        void AddOre(const MiningGeoPoly &ore, bool status = false)
        {
            mapper_.add(ore.geometry());
            if (status)
                mapper_.map(ore.geometry(), "fill-opacity:0.1;fill:#00008B;stroke:#00008B;stroke-width:0.85");
            else
                mapper_.map(ore.geometry(), "fill-opacity:0.05;fill:#51D1F6;stroke:#51D1F6;stroke-width:0.80");
        }
        void AddBlock(const MiningBlock &block, bool paint = false)
        {
            mapper_.add(block.geometry());
            if (paint)
            {
                float block_cu = block.cu();
                if (block_cu < 0.1)
                    mapper_.map(block.geometry(), "fill-opacity:0.35;fill:#C9C9C9;stroke:#C9C9C9;stroke-width:0.75");
                else if (block_cu == 0.1)
                    mapper_.map(block.geometry(), "fill-opacity:0.35;fill:#D9E1F2;stroke:#D9E1F2;stroke-width:0.75");
                else if (block_cu < 0.2)
                    mapper_.map(block.geometry(), "fill-opacity:0.35;fill:#8EA9DB;stroke:#8EA9DB;stroke-width:0.75");
                else if (block_cu == 0.2)
                    mapper_.map(block.geometry(), "fill-opacity:0.35;fill:#C6E0B4;stroke:#C6E0B4;stroke-width:0.75");
                else if (block_cu < 0.5)
                    mapper_.map(block.geometry(), "fill-opacity:0.35;fill:#00CC66;stroke:#00CC66;stroke-width:0.75");
                else if (block_cu == 0.5)
                    mapper_.map(block.geometry(), "fill-opacity:0.35;fill:#F8CBAD;stroke:#F8CBAD;stroke-width:0.75");
                else if (block_cu < 1)
                    mapper_.map(block.geometry(), "fill-opacity:0.35;fill:#C65911;stroke:#C65911;stroke-width:0.75");
                else if (block_cu == 1)
                    mapper_.map(block.geometry(), "fill-opacity:0.35;fill:#C65911;stroke:#C65911;stroke-width:0.75");
                else
                    mapper_.map(block.geometry(), "fill-opacity:0.35;fill:#8C004B;stroke:#8C004B;stroke-width:0.75");
            }
            else
                mapper_.map(block.geometry(), "fill-opacity:0.10;fill:#8C004B;stroke:#8C004B;stroke-width:0.5");
        }
        /*
        void AddSlice(const MiningSlice &slice)
        {
            mapper_.add(slice.geometry());
            mapper_.map(slice.geometry(), "fill-opacity:0;fill:#DC5300;stroke:#DC5300;stroke-width:0.5");
        }
        void AddFinger(const Finger &finger)
        {
            for (const MiningBlock &block : finger.blocks())
                AddBlock(block);
            for (const MiningGeoPoly &ore : finger.ores())
                AddOre(ore);
            for (const MiningSlice &slice : finger.slices())
                AddSlice(slice);

            AddCut(finger.cut());
        }
        void AddFingers(const std::vector<Finger> &fingers)
        {
            for (const auto &finger : fingers)
            {
                for (const MiningBlock &block : finger.blocks())
                    AddBlock(block);
            }
            for (const auto &finger : fingers)
            {
                for (const MiningSlice &slice : finger.slices())
                    AddSlice(slice);
            }
            for (const auto &finger : fingers)
            {
                for (const MiningVoid &void_ore : finger.void_ores())
                    AddVoidOre(void_ore);
            }
            for (const auto &finger : fingers)
            {
                for (const MiningGeoPoly &ore : finger.ores())
                    AddOre(ore);
            }
            for (const auto &finger : fingers)
            {
                AddCut(finger.cut(), finger.cut().extractionDay());
            }
        }
        */
    };
}