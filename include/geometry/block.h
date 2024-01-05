#pragma once

#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>

namespace bg = boost::geometry;

namespace astay
{
    class MiningBlock
    {
    private:
        typedef double float_t;
        typedef bg::model::d2::point_xy<float_t> point_xy;
        typedef bg::model::box<point_xy> box_t;
        box_t geometry_;
        int id_ = -1;
        int level_ = -1;
        float_t volume_ = -1;
        float_t density_ = -1;
        float_t cu_ = -1;
        std::string ore_type_;

    public:
        void corners(const point_xy &p_min, const point_xy &p_max)
        {
            bg::set<bg::min_corner, 0>(geometry_, p_min.x());
            bg::set<bg::min_corner, 1>(geometry_, p_min.y());
            bg::set<bg::max_corner, 0>(geometry_, p_max.x());
            bg::set<bg::max_corner, 1>(geometry_, p_max.y());
        }
        void set_id(int id)
        {
            id_ = id;
        }
        void set_level(int level)
        {
            level_ = level;
        }
        void set_volume(float_t volume)
        {
            volume_ = volume;
        }
        void set_density(float_t density)
        {
            density_ = density;
        }
        void set_cu(float_t cu)
        {
            cu_ = cu;
        }
        void set_ore_type(std::string ore_type)
        {
            ore_type_ = ore_type;
        }
        void clear()
        {
            bg::clear(geometry_);
            volume_ = -1;
            density_ = -1;
            cu_ = -1;
        }
        int id() const
        {
            return id_;
        }
        int level() const
        {
            return level_;
        }
        float_t volume() const
        {
            return volume_;
        }
        float_t density() const
        {
            return density_;
        }
        float_t cu() const
        {
            return cu_;
        }
        std::string ore_type() const
        {
            return ore_type_;
        }
        const box_t &geometry() const
        {
            return geometry_;
        }
        void correct()
        {
            bg::correct(geometry_);
        }

        friend std::ostream &operator<<(std::ostream &o, const MiningBlock &block)
        {
            return o << "[" << block.id_ << "] | " << bg::dsv(block.geometry_)
                     << " | vol:" << block.volume_
                     << ", cu:" << block.cu_
                     << ", density:" << block.density_;
        }
    };
}