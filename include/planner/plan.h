#pragma once

#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>

namespace bg = boost::geometry;

namespace astay
{
class MiningPlan
{
private:
    typedef double float_t;
    float_t xOffset_, yOffset_;
    std::vector<float_t> x_, y_;
    std::vector<float_t> levels_;
    std::vector<std::string> loaders_;
    int size_ = 0;

public:
    void add_row(float_t x, float_t y, float_t level, std::string &loader)
    {
        x_.push_back(x);
        y_.push_back(y);
        levels_.push_back(level);
        loaders_.push_back(loader);
        size_++;
    }
    void add_offset(float_t x_offset, float_t y_offset)
    {
        xOffset_ = x_offset;
        yOffset_ = y_offset;
    }
    const std::vector<float_t> &x() const
    {
        return x_;
    }
    const std::vector<float_t> &y() const
    {
        return y_;
    }
    const std::vector<float_t> &levels() const
    {
        return levels_;
    }
    const std::vector<std::string> &loaders() const
    {
        return loaders_;
    }
    int size()
    {
        return size_;
    }
    void clear()
    {
        size_ = 0;
        x_.clear();
        y_.clear();
        levels_.clear();
        loaders_.clear();
    }
    friend std::ostream &operator<<(std::ostream &o, const MiningPlan &plan)
    {
        int idx = 0;
        for (auto &loader : plan.loaders())
        {
            o << "[" << loader << "] | ("
              << plan.x()[idx] << ", "
              << plan.y()[idx] << ", "
              << plan.levels()[idx] << ")" << std::endl;
              idx++;
        }
        return o;
    }
};
}