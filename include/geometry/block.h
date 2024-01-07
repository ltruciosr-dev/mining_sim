#pragma once

#include "geometry/polygon.h"

namespace bg = boost::geometry;

namespace astay
{
    class MiningBlock : public Polygon
    {
    private:
        int type_{-1};
        float_t topo_{-1}, density_{-1};
        float_t au_{-1}, ag_{-1}, cu_{-1};
        float_t tcm_{-1}, s_{-1};

    public:
        void setType(int type)
        {
            type_ = type;
        }
        void setTopo(float_t topo)
        {
            topo_ = topo;
        }
        void setDensity(float_t density)
        {
            density_ = density;
        }
        void setAu(float_t au)
        {
            au_ = au;
        }
        void setAg(float_t ag)
        {
            ag_ = ag;
        }
        void setCu(float_t cu)
        {
            cu_ = cu;
        }
        void setTcm(float_t tcm)
        {
            tcm_ = tcm;
        }
        void setS(float_t s)
        {
            s_ = s;
        }
        int type() const
        {
            return type_;
        }
        float_t topo() const
        {
            return topo_;
        }
        float_t density() const
        {
            return density_;
        }
        float_t au() const
        {
            return au_;
        }
        float_t ag() const
        {
            return ag_;
        }
        float_t cu() const
        {
            return cu_;
        }
        float_t tcm() const
        {
            return tcm_;
        }
        float_t s() const
        {
            return s_;
        }
        friend std::ostream &operator<<(std::ostream &o, const MiningBlock &block)
        {
            return o << "[" << block.id_ << "] | " << bg::dsv(block.geometry_)
                     << " | TOPO%:" << block.topo_
                     << ", cu:" << block.cu_
                     << ", density:" << block.density_;
        }
    };

    class SliderBlock : public Polygon
    {
    private:
    public:
    
    }
}