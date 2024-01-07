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
        int level_, id_;

    public:
        void add(const point_xy &p)
        {
            bg::append(geometry_, p);
        }
        void setLevel(int level)
        {
            level_ = level;
        }
        void setId(int id)
        {
            id_ = id;
        }
        void computeBox()
        {
            float_t x_min{DBL_MAX}, y_min{DBL_MAX};
            float_t x_max{-DBL_MAX}, y_max{-DBL_MAX};
        
            // Iterate through vertex to find the extreme coordinates
            for (auto it = boost::begin(bg::exterior_ring(geometry_)); it != boost::end(bg::exterior_ring(geometry_)); ++it)
            {
                float_t x = bg::get<0>(*it);
                float_t y = bg::get<1>(*it);
                // Update the overall min and max coordinates
                x_min = std::min(x, x_min);
                y_min = std::min(y, y_min);
                x_max = std::max(x, x_max);
                y_max = std::max(y, y_max);
            }
            bg::set<bg::min_corner, 0>(box_, x_min);
            bg::set<bg::min_corner, 1>(box_, y_min);
            bg::set<bg::max_corner, 0>(box_, x_max);
            bg::set<bg::max_corner, 1>(box_, y_max);
        }
        void correct()
        {
            bg::correct(geometry_);
        }
        void clear()
        {
            bg::clear(geometry_);
            bg::clear(box_);
        }
        const box_t &box() const
        {
            return box_;
        }
        const polygon_t &geometry() const
        {
            return geometry_;
        }
        int level() const
        {
            return level_;
        }
        int id() const
        {
            return id_;
        }
        int size() const
        {
            return bg::num_points(geometry_);
        }
        bool containsPolygon(const polygon_t &geometry)
        {
            if (bg::intersects(geometry_, geometry))
            {
                return true;
            }
            return false;
        }
    };

    class MiningCut : public Polygon
    {
    private:
        std::string name_;
        int extraction_day_;

    public:
        void setName(std::string &name)
        {
            name_ = name;
        }
        void setExtractionDay(int extraction_day)
        {
            extraction_day_ = extraction_day;
        }
        std::string name() const
        {
            return name_;
        }
        int extractionDay() const
        {
            return extraction_day_;
        }
        friend std::ostream &operator<<(std::ostream &o, const MiningCut &poly)
        {
            return o << "Mining Cut - ID (" << poly.id() << ") :\tVertices (" << bg::num_points(poly.geometry()) 
                     << ")\t| Level (" << poly.level() << ") | Day (" << poly.extractionDay() << ")";
        }
    };

    class MiningGeoPoly : public Polygon
    {
    private:
        std::string name_ = "NN";
        bool status_ = false;

    public:
        void setName(std::string &name)
        {
            name_ = name;
        }
        void setStatus(bool status)
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
        friend std::ostream &operator<<(std::ostream &o, const MiningGeoPoly &ore)
        {
            return o << "[" << ore.name_ << "] | " << bg::dsv(ore.box_) << " | level(" << ore.level_ << ")";
        }
    };
}