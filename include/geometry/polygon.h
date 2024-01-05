#pragma once

#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include "geometry/block.h"

namespace bg = boost::geometry;

namespace astay
{
    class Polygon
    {
    protected:
        typedef double float_t;
        typedef bg::model::d2::point_xy<float_t> point_xy;
        typedef bg::model::polygon<point_xy> polygon_t;
        typedef bg::model::box<point_xy> box_t;
        polygon_t geometry_;
        box_t box_;
        int level_;

    public:
        void add(const point_xy &p)
        {
            bg::append(geometry_, p);
        }
        void set_level(int level)
        {
            level_ = level;
        }
        void compute_box()
        {
            // Initialize box vertices
            float_t min_x(FLT_MAX), min_y(FLT_MAX), max_x(-FLT_MAX), max_y(-FLT_MAX);

            for (auto it = boost::begin(bg::exterior_ring(geometry_)); it != boost::end(bg::exterior_ring(geometry_)); ++it)
            {
                float_t x = bg::get<0>(*it);
                float_t y = bg::get<1>(*it);
                // define x
                if (x < min_x)
                    min_x = x;
                else if (x > max_x)
                    max_x = x;
                // define y
                if (y < min_y)
                    min_y = y;
                else if (y > max_y)
                    max_y = y;
            }
            bg::set<bg::min_corner, 0>(box_, min_x);
            bg::set<bg::min_corner, 1>(box_, min_y);
            bg::set<bg::max_corner, 0>(box_, max_x);
            bg::set<bg::max_corner, 1>(box_, max_y);
        }
        const box_t &box() const
        {
            return box_;
        }
        const polygon_t &geometry() const
        {
            return geometry_;
        }
        const int level() const
        {
            return level_;
        }
        void correct()
        {
            bg::correct(geometry_);
        }
        bool ContainsPolygon(const polygon_t &geometry)
        {
            if (bg::intersects(geometry_, geometry))
            {
                return true;
            }
            return false;
        }
        bool ContainsBlock(const MiningBlock &block)
        {
            if (bg::intersects(geometry_, block.geometry()))
            {
                return true;
            }
            return false;
        }
        int size() const
        {
            return bg::num_points(geometry_);
        }
        void clear()
        {
            bg::clear(geometry_);
            bg::clear(box_);
        }
        friend std::ostream &operator<<(std::ostream &o, const Polygon &poly)
        {
            return o << bg::dsv(poly.box_) << " | level(" << poly.level_ << ")";
        }
        static box_t GetBoxFromPolygon(const polygon_t &polygon)
        {
            box_t box;
            float_t min_x(FLT_MAX), min_y(FLT_MAX), max_x(-FLT_MAX), max_y(-FLT_MAX);
            for (auto it = boost::begin(bg::exterior_ring(polygon)); it != boost::end(bg::exterior_ring(polygon)); ++it)
            {
                float_t x = bg::get<0>(*it);
                float_t y = bg::get<1>(*it);
                // define x
                if (x < min_x)
                    min_x = x;
                else if (x > max_x)
                    max_x = x;
                // define y
                if (y < min_y)
                    min_y = y;
                else if (y > max_y)
                    max_y = y;
            }
            bg::set<bg::min_corner, 0>(box, min_x);
            bg::set<bg::min_corner, 1>(box, min_y);
            bg::set<bg::max_corner, 0>(box, max_x);
            bg::set<bg::max_corner, 1>(box, max_y);
            return box;
        }
    };

    class MiningCut : public Polygon
    {
    private:
        std::string name_;
        int extraction_day_;
        int id_;

    public:
        void set_name(std::string &name)
        {
            name_ = name;
        }
        void set_extraction_day(int extraction_day)
        {
            extraction_day_ = extraction_day;
        }
        void set_id(int id)
        {
            id_ = id;
        }
        std::string name() const
        {
            return name_;
        }
        int extraction_day() const
        {
            return extraction_day_;
        }
        int id() const
        {
            return id_;
        }
        friend std::ostream &operator<<(std::ostream &o, const MiningCut &poly)
        {
            return o << "vertices (" << bg::num_points(poly.geometry_) << ") | level (" << poly.level_ << ")";
        }
    };

    class MiningOre : public Polygon
    {
    private:
        std::string name_ = "NN";
        bool status_ = false;

    public:
        void set_name(std::string &name)
        {
            name_ = name;
        }
        void set_status(bool status)
        {
            status_ = status;
        }
        std::string name() const
        {
            return name_;
        }
        bool status() const
        {
            return status_;
        }
        void clear()
        {
            bg::clear(geometry_);
            bg::clear(box_);
            name_ = "NN";
            status_ = false;
        }
        friend std::ostream &operator<<(std::ostream &o, const MiningOre &ore)
        {
            return o << "[" << ore.name_ << "] | " << bg::dsv(ore.box_) << " | level(" << ore.level_ << ")";
        }
    };

    class MiningVoid : public Polygon
    {
    public:
        void set_geometry(polygon_t geometry, float_t level)
        {
            geometry_ = geometry;
            compute_box();
            set_level(level);
            correct();
        }
        friend std::ostream &operator<<(std::ostream &o, const MiningVoid &ore)
        {
            return o << "[void_cut] | " << bg::dsv(ore.box_) << " | level(" << ore.level_ << ")";
        }
    };

    class MiningSlice : public Polygon
    {
    public:
    };
}