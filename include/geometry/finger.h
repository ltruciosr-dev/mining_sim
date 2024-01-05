#pragma once

#include <vector>
#include <boost/geometry.hpp>
#include <boost/foreach.hpp>
#include "geometry/polygon.h"
#include "geometry/block.h"

namespace bg = boost::geometry;

namespace astay
{
    class Finger
    {
    private:
        typedef double float_t;
        typedef bg::model::d2::point_xy<float_t> point_xy;
        typedef bg::model::polygon<point_xy> polygon_t;
        MiningCut cut_;
        std::vector<MiningBlock> blocks_;
        std::vector<MiningOre> ores_;
        std::vector<MiningSlice> slices_;
        std::vector<MiningVoid> void_ores_;
        int extraction_day_;
        int id_;

    public:
        void set_cut(MiningCut cut)
        {
            cut_ = cut;
        }
        void set_day(int extraction_day)
        {
            extraction_day_ = extraction_day;
        }
        void set_id(int id)
        {
            id_ = id;
        }
        void AddBlocks(std::vector<MiningBlock> blocks)
        {
            blocks_ = blocks;
        }
        void AddOre(MiningOre ore)
        {
            ores_.push_back(ore);
        }
        void AddSlice(MiningSlice slice)
        {
            slices_.push_back(slice);
        }
        bool ContainsOre(MiningOre ore)
        {
            if (bg::intersects(cut_.box(), ore.box()))
            {
                if (bg::intersects(cut_.geometry(), ore.geometry()))
                    return true;
            }
            return false;
        }
        bool ContainsBlock(MiningBlock block)
        {
            if (bg::intersects(cut_.geometry(), block.geometry()))
            {
                return true;
            }
            return false;
        }
        bool ContainsSlice(MiningSlice slice)
        {
            if (bg::intersects(cut_.box(), slice.box()))
            {
                if (bg::intersects(cut_.geometry(), slice.geometry()))
                    return true;
            }
            return false;
        }
        int day() const
        {
            return extraction_day_;
        }
        int id() const
        {
            return id_;
        }
        const MiningCut &cut() const
        {
            return cut_;
        }
        const std::vector<MiningBlock> &blocks() const
        {
            return blocks_;
        }
        const std::vector<MiningOre> &ores() const
        {
            return ores_;
        }
        const std::vector<MiningSlice> &slices() const
        {
            return slices_;
        }
        const std::vector<MiningVoid> &void_ores() const
        {
            return void_ores_;
        }
        void ComputeSlices(float_t slice_dist, point_xy dir_vec)
        {
            // Get limits by given direction
            point_xy ref_point{0, 0};
            auto limits = ComputeLimitsByDirection(dir_vec, ref_point);
            float_t u_min{limits[0]}, u_max{limits[1]};
            float_t v_min{limits[2]}, v_max{limits[3]};

            // Generate slices
            int n_slices = ceil((u_max - u_min) / slice_dist);
            point_xy orth_vec{-bg::get<1>(dir_vec), bg::get<0>(dir_vec)}; // Orthogonal vector
            point_xy lslices_base{bg::get<0>(ref_point) + u_min * bg::get<0>(dir_vec),
                                  bg::get<1>(ref_point) + u_min * bg::get<1>(dir_vec)};
            point_xy lslices_lower{bg::get<0>(lslices_base) + v_min * bg::get<0>(orth_vec),
                                   bg::get<1>(lslices_base) + v_min * bg::get<1>(orth_vec)};
            point_xy lslices_upper{bg::get<0>(lslices_base) + v_max * bg::get<0>(orth_vec),
                                   bg::get<1>(lslices_base) + v_max * bg::get<1>(orth_vec)};

            for (int i = 0; i < n_slices; i++)
            {
                MiningSlice slice;
                point_xy rslices_upper{bg::get<0>(lslices_upper) + slice_dist * bg::get<0>(dir_vec),
                                       bg::get<1>(lslices_upper) + slice_dist * bg::get<1>(dir_vec)};
                point_xy rslices_lower{bg::get<0>(lslices_lower) + slice_dist * bg::get<0>(dir_vec),
                                       bg::get<1>(lslices_lower) + slice_dist * bg::get<1>(dir_vec)};
                slice.add(lslices_upper);
                slice.add(lslices_lower);
                slice.add(rslices_lower);
                slice.add(rslices_upper);
                slice.add(lslices_upper);
                slice.compute_box();
                slice.correct();
                slices_.push_back(slice);
                // Update
                lslices_upper = rslices_upper;
                lslices_lower = rslices_lower;
            }
        }
        std::array<float_t, 4> ComputeLimitsByDirection(point_xy &dir_vec, point_xy &ref_point) const
        {
            float_t u_min{FLT_MAX}, u_max{-FLT_MAX}; // parallel to dir_vec
            float_t v_min{FLT_MAX}, v_max{-FLT_MAX}; // perpendicular to dir_vec
            // Init variables
            float_t u_ref = bg::get<0>(ref_point);
            float_t v_ref = bg::get<1>(ref_point);
            point_xy orth_vec{-bg::get<1>(dir_vec), bg::get<0>(dir_vec)};
            
            // Identify the nearest and the farthest vertices with the provided direction.
            for (auto it = boost::begin(bg::exterior_ring(cut_.geometry())); it != boost::end(bg::exterior_ring(cut_.geometry())); ++it) // Iterator of polygon vertices
            {
                float_t u = bg::get<0>(*it);
                float_t v = bg::get<1>(*it);
                point_xy vector{u - u_ref, v - v_ref};
                float_t u_dist = bg::dot_product(dir_vec, vector);
                float_t v_dist = bg::dot_product(orth_vec, vector);
                // Get 'u' limits
                if (u_dist < u_min)
                {
                    u_min = u_dist;
                }
                else if (u_dist > u_max)
                {
                    u_max = u_dist;
                }
                // Get 'v' limits
                if (v_dist < v_min)
                {
                    v_min = v_dist;
                }
                else if (v_dist > v_max)
                {
                    v_max = v_dist;
                }
            }
            std::array<float_t, 4> limits = {u_min, u_max, v_min, v_max};
            return limits;
        }
        void FixMissingOre(float_t min_area)
        {
            std::vector<polygon_t> cut_diff;
            cut_diff.push_back(cut_.geometry());
            for (const auto &ore : ores_)
            {
                std::vector<polygon_t> new_cut_diff;
                for (auto &diff : cut_diff)
                {
                    bg::correct(diff);
                    boost::geometry::difference(diff, ore.geometry(), new_cut_diff);
                }
                cut_diff = new_cut_diff;
            }
            float_t level = cut_.level();
            void_ores_.clear();
            for (auto &diff : cut_diff)
            {
                if (bg::area(diff) > min_area)
                {
                    MiningVoid missing_ore;
                    missing_ore.set_geometry(diff, level);
                    void_ores_.push_back(missing_ore);
                }
            }
        }
        void clear()
        {
            cut_.clear();
            blocks_.clear();
            ores_.clear();
            slices_.clear();
            void_ores_.clear();
            extraction_day_ = -1;
        }
        friend std::ostream &operator<<(std::ostream &o, const Finger &finger)
        {
            int idx = 1;
            o << "- Polygon: (" << finger.cut_.size() << " vertices)" << '\n'
              << "- Blocks: (" << finger.blocks_.size() << " units)" << '\n';
            o << "- Ores: (" << finger.ores_.size() << " units)" << '\n';
            o << "- Missing Ores (" << finger.void_ores_.size() << " units)" << '\n';
            for (const auto &void_ore : finger.void_ores_)
            {
                o << "  > missing ore " << idx++ << ": " << bg::area(void_ore.geometry()) << '\n';
            }
            o << "- Slices: (" << finger.slices_.size() << " units)" << '\n';
            return o;
        }
    };
}