#pragma once

#include <math.h>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include "SQLiteCpp/SQLiteCpp.h"
#include "geometry/polygon.h"
#include "geometry/block.h"
#include "geometry/finger.h"

namespace astay
{
    class ModelBlock
    {
    private:
        typedef double float_t;
        typedef bg::model::d2::point_xy<float_t> point_xy;
        SQLite::Database db_;
        float_t num_cells_x_, num_cells_y_, num_cells_z_;
        float_t orig_x_, orig_y_, orig_z_;
        float_t block_size_;
        float_t offset_x_{0}, offset_y_{0}, offset_z_{0};

    public:
        ModelBlock(const std::string &db_name) : db_(db_name)
        {
        }
        void set_offset(float_t offset_x = 0, float_t offset_y = 0, float_t offset_z = 0)
        {
            offset_x_ = offset_x;
            offset_y_ = offset_y;
            offset_z_ = offset_z;
        }
        void read_setup()
        {
            SQLite::Statement query(db_, "SELECT * from ProjectModel");
            while (query.executeStep())
            {
                num_cells_x_ = query.getColumn(3).getDouble();
                num_cells_y_ = query.getColumn(4).getDouble();
                num_cells_z_ = query.getColumn(5).getDouble();
                orig_x_ = query.getColumn(6).getDouble() - offset_x_;
                orig_y_ = query.getColumn(7).getDouble() - offset_y_;
                orig_z_ = query.getColumn(8).getDouble() - offset_z_;
                block_size_ = query.getColumn(9).getDouble();
            }
        }
        std::vector<MiningBlock> get_blocks(MiningCut &cut)
        {
            std::vector<MiningBlock> blocks;
            float_t poly_xmin = bg::get<bg::min_corner, 0>(cut.box());
            float_t poly_ymin = bg::get<bg::min_corner, 1>(cut.box());
            float_t poly_xmax = bg::get<bg::max_corner, 0>(cut.box());
            float_t poly_ymax = bg::get<bg::max_corner, 1>(cut.box());
            float_t poly_level = cut.level();

            // Compute limits
            int delta_xmin = floor((poly_xmin - orig_x_) / block_size_);
            int delta_ymin = floor((poly_ymin - orig_y_) / block_size_);
            int delta_xmax = floor((poly_xmax - orig_x_) / block_size_);
            int delta_ymax = floor((poly_ymax - orig_y_) / block_size_);
            int delta_level = floor((poly_level - orig_z_) / block_size_);

            // Compute blocks
            MiningBlock block;
            int id_x, id_y, id_z;
            int id; // sum(id_x, id_y, id_z)
            std::string query_str = "SELECT Volumen, Density, Cu, OreSort FROM ModelBlock WHERE rowid = ";

            id_z = delta_level * (num_cells_x_ * num_cells_y_);
            for (int delta_y = delta_ymin; delta_y <= delta_ymax; delta_y++)
            {
                id_y = delta_y * num_cells_x_;
                for (int delta_x = delta_xmin; delta_x <= delta_xmax; delta_x++)
                {
                    id_x = delta_x + 1;
                    id = id_x + id_y + id_z; // ID to read SQL database
                    float_t block_xmin = orig_x_ + delta_x * block_size_;
                    float_t block_ymin = orig_y_ + delta_y * block_size_;
                    point_xy p_min{block_xmin, block_ymin};
                    point_xy p_max{block_xmin + block_size_, block_ymin + block_size_};
                    block.corners(p_min, p_max);
                    if (cut.ContainsBlock(block))
                    {
                        block.set_id(id);
                        block.correct();

                        // SQLite Interaction
                        std::string ore_type;
                        float_t volume, density, cu;
                        SQLite::Statement query(db_, query_str + std::to_string(id));
                        while (query.executeStep())
                        {
                            volume = query.getColumn(0).getDouble();
                            density = query.getColumn(1).getDouble();
                            cu = query.getColumn(2).getDouble();
                            ore_type = query.getColumn(3).getString();
                        }
                        block.set_level(poly_level);
                        block.set_density(density);
                        block.set_volume(volume);
                        block.set_cu(cu);
                        block.set_ore_type(ore_type);
                        blocks.push_back(block);
                    }
                    block.clear();
                }
            }
            return blocks;
        }
        std::vector<MiningBlock> get_blocks(Finger &finger)
        {
            auto cut = finger.cut();
            return get_blocks(cut);
        }
        friend std::ostream &operator<<(std::ostream &o, const ModelBlock &model_block)
        {
            o << "- Origen: ";
            o << model_block.orig_x_ << ", "
              << model_block.orig_y_ << ", "
              << model_block.orig_z_ << std::endl;
            o << "- Offset: ";
            o << model_block.offset_x_ << ", "
              << model_block.offset_y_ << ", "
              << model_block.offset_z_ << ", " << std::endl;
            o << "- Size Block: ";
            o << model_block.block_size_ << std::endl;
            o << "- # Cells: (";
            o << model_block.num_cells_x_ << ", "
              << model_block.num_cells_y_ << ", "
              << model_block.num_cells_z_ << ")" << std::endl;
            return o;
        }
    };
}