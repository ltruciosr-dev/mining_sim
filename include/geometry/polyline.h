#pragma once
#include <iostream>
#include <fstream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

namespace bg = boost::geometry;

namespace astay
{
    class Polyline
    {
    protected:
        typedef double float_t;
        typedef bg::model::d2::point_xy<float_t> point_t;
        typedef bg::model::linestring<point_t> linestring_t;
        linestring_t linestring_;
        point_t centroid_;
        int level_{-1}, id_{-1};        
    public:
        void add(const point_t &p)
        {
            bg::append(linestring_, p);
        }
        void setLevel(int level)
        {
            level_ = level;
        }
        void setId(int id)
        {
            id_ = id;
        }
        void computeCentroid()
        {
            bg::centroid(linestring_, centroid_);
        }
        void clear()
        {
            bg::clear(linestring_);
        }
        const linestring_t &linestring() const
        {
            return linestring_;
        }
        int level() const
        {
            return level_;
        }
        int id() const
        {
            return id_;
        }
        point_t centroid() const
        {
            return centroid_;
        }
        int size() const
        {
            return bg::length(linestring_);
        }
        friend std::ostream &operator<<(std::ostream &o, const Polyline &polyline)
        {
            return o << "Polyline - ID (" << polyline.id() << ") : " << polyline.size() << " pts. \t| level (" << polyline.level_ << ")";
        }
    };
}