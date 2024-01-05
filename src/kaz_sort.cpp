#include <numeric>
#include <algorithm>
#include "rapidcsv.h"
#include "kaz_sort.h"
#include "geometry/polygon.h"
#include "geometry/drawer.h"

#include "SQLiteCpp/Database.h"
#include "SQLiteCpp/Transaction.h"

namespace astay
{
    KazSort::KazSort(std::string db_name) : db_(db_name), transaction_(db_)
    {
    }

    // -- set
    void KazSort::setProjectID(int project_id)
    {
        project_id_ = project_id;
    }

    void KazSort::setProjectDB(const std::string &db_name)
    {
        project_db_ = db_name;
    }

    void KazSort::setModelBlocksDB(const std::string &db_name)
    {
        modelblocks_db_ = db_name;
    }

    void KazSort::setSVGDir(const std::string &svg_dir)
    {
        svg_dir_ = svg_dir;
    }

    // -- read
    void KazSort::ReadCuts(bool verbose)
    {
        SQLite::Database db(project_db_, SQLite::OPEN_READONLY);
        std::vector<float_t> x_cut, y_cut, z_cut;
        std::vector<int> vert_num, extraction_day;

        std::string query_cuts = "SELECT rowid,X,Y,Z,Vertices,LayerID from MiningCuts WHERE Project_ID=" + std::to_string(project_id_) + " ORDER BY IndexPoly ASC, Vertices ASC";
        SQLite::Statement query(db, query_cuts);
        while (query.executeStep())
        {
            int rowid = query.getColumn(0).getInt();
            float_t x = query.getColumn(1).getDouble();
            float_t y = query.getColumn(2).getDouble();
            float_t z = query.getColumn(3).getDouble();
            int vertice = query.getColumn(4).getInt();
            int day = query.getColumn(5).getInt();
            // Populate vectors
            row_ids_.push_back(rowid);
            x_cut.push_back(x);
            y_cut.push_back(y);
            z_cut.push_back(z);
            vert_num.push_back(vertice);
            extraction_day.push_back(day);
        }

        // validate
        if (row_ids_.size() < 3)
        {
            if (verbose)
                std::cout << "[ERROR] ReadCuts | Doesn't have any lines." << std::endl;
            return;
        }
        if (verbose)
            std::cout << "[DEBUG] ReadCuts | Has " << vert_num.size() << " lines." << std::endl;

        // Offset
        x_offset_ = std::accumulate(x_cut.begin(), x_cut.end(), 0) / x_cut.size();
        y_offset_ = std::accumulate(y_cut.begin(), y_cut.end(), 0) / y_cut.size();
        for (float_t &x : x_cut)
            x -= x_offset_;
        for (float_t &y : y_cut)
            y -= y_offset_;

        // Extract
        int p_vert_n = vert_num[0] - 1;
        MiningCut cut;
        for (size_t i = 0; i < vert_num.size(); i++)
        {
            int vert_n = vert_num[i];
            if (vert_n != p_vert_n + 1) // if true, reading new cut
            {
                cut.compute_box();
                cut.set_level(z_cut[i - 1]);
                cut.set_extraction_day(extraction_day[i - 1]);
                cut.correct();
                cuts_.push_back(cut);
                cut.clear();
            }
            point_xy p{x_cut[i], y_cut[i]};
            cut.add(p);
            p_vert_n = vert_n;
        }
        cut.compute_box();
        cut.set_level(z_cut.back());
        cut.set_extraction_day(extraction_day.back());
        cut.correct();
        cuts_.push_back(cut);

        if (verbose)
        {
            std::cout << "Offset: (" << x_offset_ << ", " << y_offset_ << ")" << std::endl;
            std::cout << "Number of Cuts: " << cuts_.size() << std::endl;
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
            for (auto cut : cuts_)
            {
                std::cout << "- " << cut << std::endl;
            }
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
        }
    }

    void KazSort::ReadOres(bool verbose)
    {
        // store the cut's index the belogs to each ore_level
        std::map<int, std::vector<int>> ore_setup; // (ore_level) -> [cut_idx_0, cut_idx_1, ...]
        int idx = 0;
        for (MiningCut &cut : cuts_)
            ore_setup[cut.level()].push_back(idx++);

        // extract ore_polygons that are retrieved in ore_setup
        SQLite::Database db(modelblocks_db_, SQLite::OPEN_READONLY);
        for (auto const &set : ore_setup)
        {
            // extract all the ore polygons that belongs to each level
            int level = set.first;

            std::vector<float_t> x_ore, y_ore, z_ore;
            std::vector<int> vert_num;
            std::vector<std::string> ore_name;

            std::string query_ore = "SELECT X,Y,Z,Number,Element_Name from OrePoly WHERE Z=" + std::to_string(level) + " ORDER BY Element_Name ASC, Number ASC";
            SQLite::Statement query(db, query_ore);
            while (query.executeStep())
            {
                float_t x = query.getColumn(0).getDouble();
                float_t y = query.getColumn(1).getDouble();
                float_t z = query.getColumn(2).getDouble();
                int vertice = query.getColumn(3).getInt();
                std::string name = query.getColumn(4).getString();
                // Populate vectors
                x_ore.push_back(x);
                y_ore.push_back(y);
                z_ore.push_back(z);
                vert_num.push_back(vertice);
                ore_name.push_back(name);
            }

            // validate
            if (vert_num.size() < 3)
            {
                if (verbose)
                    std::cout << "[ERROR] ReadOres | level " << level << " doesn't have any lines." << std::endl;
                continue;
            }
            if (verbose)
                std::cout << "[DEBUG] ReadOres | level " << level << " have " << vert_num.size() << " lines." << std::endl;

            // substract offset
            for (float_t &x : x_ore)
                x -= x_offset_;
            for (float_t &y : y_ore)
                y -= y_offset_;

            // init extraction
            int p_vert_n = vert_num[0] - 1;
            MiningOre ore;
            for (size_t i = 0; i < vert_num.size(); i++)
            {
                int vert_n = vert_num[i];
                if (vert_n != p_vert_n + 1) // if true, reading new ore
                {
                    ore.compute_box();
                    ore.set_name(ore_name[i - 1]);
                    ore.set_level(z_ore[i - 1]);
                    ore.correct();
                    for (int idx : set.second)
                    {
                        if (cuts_[idx].ContainsPolygon(ore.geometry()))
                        {
                            ore.set_status(true); // false by default
                            break;
                        }
                    }
                    ores_by_level_[ore.level()].push_back(ore);
                    ore.clear();
                }
                point_xy p{x_ore[i], y_ore[i]};
                ore.add(p);
                p_vert_n = vert_n;
            }
            ore.compute_box();
            ore.set_name(ore_name.back());
            ore.set_level(z_ore.back());
            ore.correct();
            for (int idx : set.second)
            {
                if (cuts_[idx].ContainsPolygon(ore.geometry()))
                {
                    ore.set_status(true);
                    break;
                }
            }
            ores_by_level_[ore.level()].push_back(ore);
        }

        if (verbose)
        {
            std::cout << "Offset: (" << x_offset_ << ", " << y_offset_ << ")" << std::endl;
            std::cout << "Ore has (" << ores_by_level_.size() << ") levels." << std::endl;
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
            for (const auto &ores_level : ores_by_level_)
            {
                std::cout << "- level [" << ores_level.first << "] | polygons [" << ores_level.second.size() << "]" << std::endl;
            }
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
        }
    }

    void KazSort::ReadBlocks(bool verbose)
    {
        ModelBlock model_block_(modelblocks_db_);
        model_block_.set_offset(x_offset_, y_offset_);
        model_block_.read_setup();
        std::set<int> blocks_id;
        for (auto &cut : cuts_)
        {
            std::vector<MiningBlock> temp_blocks = model_block_.get_blocks(cut);
            for (MiningBlock block : temp_blocks)
            {
                int id = block.id();
                if (blocks_id.find(id) != blocks_id.end())
                    continue;
                else
                {
                    blocks_id.insert(id);
                    blocks_by_level_[block.level()].push_back(block);
                }
            }
        }
    }

    void KazSort::ReadPlan(bool verbose)
    {
        // Clear
        plan_.clear();
        // Extract from DB
        SQLite::Database db(project_db_, SQLite::OPEN_READONLY);
        std::string query_plan = "SELECT Source_,X_Source,Y_Source,Z_Source from Vista_Material WHERE Project_Id=" + std::to_string(project_id_);
        SQLite::Statement query(db, query_plan);
        while (query.executeStep())
        {
            std::string source = query.getColumn(0).getString();
            std::string loader = source.substr(0, source.find('_'));
            float_t x = query.getColumn(1).getDouble() - x_offset_;
            float_t y = query.getColumn(2).getDouble() - y_offset_;
            float_t z = query.getColumn(3).getDouble();
            // Populate vectors
            plan_.add_row(x, y, z, loader);
        }
        if (verbose)
        {
            std::cout << "Number of rows: " << plan_.size() << std::endl;
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
            std::cout << plan_ << std::endl;
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
        }
    }

    // -- compute
    void KazSort::SortCuts(bool verbose)
    {
        // Clear vectors
        cuts_by_loader_.clear();
        row_loaders_.clear();
        row_cuts_.clear();
        for (auto &cut : cuts_)
        {
            // Extract cut centroid
            point_xy centroid;
            bg::centroid(cut.geometry(), centroid);

            int cut_vert = bg::num_points(cut.geometry());

            // Get the loader whose some position is nearest to the cut
            std::string loader;
            float_t min_dist = FLT_MAX;
            for (int i = 0; i < plan_.size(); i++)
            {
                // Distance between cut centroid and each plan
                float dist = sqrt(pow(plan_.x()[i] - centroid.x(), 2) +
                                  pow(plan_.y()[i] - centroid.y(), 2));
                // Store neariest plan loader
                if (dist < min_dist)
                {
                    min_dist = dist;
                    loader = plan_.loaders()[i];
                }
            }
            int cut_idx = cuts_by_loader_[loader].size() + 1;
            cut.set_id(cut_idx);
            // Populate loader per vertice
            row_loaders_.insert(row_loaders_.end(), cut_vert, loader);
            row_cuts_.insert(row_cuts_.end(), cut_vert, cut_idx);
            // Populate a dict of cuts by loader
            cuts_by_loader_[loader].push_back(cut);
        }
        if (verbose)
        {
            std::cout << "Number of Loaders: " << cuts_by_loader_.size() << std::endl;
            std::cout << "Number of Cuts: " << cuts_.size() << std::endl;
            std::cout << "Number of Vertices: " << row_loaders_.size() << std::endl;
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
            for (auto const &loader_cuts : cuts_by_loader_)
            {
                std::cout << loader_cuts.first << ":" << std::endl;
                for (auto const &cut : loader_cuts.second)
                {
                    point_xy cut_centroid;
                    bg::centroid(cut.geometry(), cut_centroid);
                    std::cout << "- " << cut
                              << " | centroid ("
                              << cut_centroid.x() << ","
                              << cut_centroid.y() << ")" << std::endl;
                }
            }
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
        }
    }

    void KazSort::ValidateDB(bool verbose)
    {
        SQLite::Database db(project_db_, SQLite::OPEN_READONLY);
        std::string query_direction = "SELECT COUNT(Project_ID) FROM MiningDirection WHERE Project_ID=" + std::to_string(project_id_);
        SQLite::Statement query(db, query_direction);
        int count;
        while (query.executeStep())
        {
            count = query.getColumn(0).getInt();
        }
        if (count > 0)
            is_sorted_ = true;
        else
            is_sorted_ = false;

        if (verbose)
            std::cout << "Validation Count: " << count << std::endl;
    }

    void KazSort::PopulateDB(bool verbose)
    {
        if (is_sorted_)
        {
            if (verbose)
                std::cout << "Database already updated" << std::endl;
        }
        else
        {
            SQLite::Database db(project_db_, SQLite::OPEN_READWRITE);
            SQLite::Transaction transaction(db);
            // Populate MiningCuts
            int idx = 0;
            for (const auto &loader : row_loaders_)
            {
                int rowid = row_ids_[idx];
                int cutid = row_cuts_[idx];
                // populate 'Loader' and 'Cut_ID'
                std::string query_loader = "UPDATE MiningCuts SET Loader='" + loader +
                                           "' WHERE rowid=" + std::to_string(rowid);
                std::string query_cut = "UPDATE MiningCuts SET Cut_ID='" + std::to_string(cutid) +
                                        "' WHERE rowid=" + std::to_string(rowid);
                db.exec(query_loader);
                db.exec(query_cut);
                idx++;
                if (verbose)
                {
                    std::cout << "Loader query:" << query_loader << std::endl;
                    std::cout << "Cut query:" << query_cut << std::endl;
                }
            }
            // Populate MiningDirection
            for (const auto &loader_cuts : cuts_by_loader_)
            {
                std::string shovel = loader_cuts.first;
                for (const auto &cut : loader_cuts.second)
                {
                    std::string query_direction = "INSERT INTO MiningDirection (Shovel,Cut_ID,Dir_X,Dir_Y,Project_ID) VALUES ('" +
                                                  shovel + "'," + std::to_string(cut.id()) + ",0.0,0.0," + std::to_string(project_id_) + ")";
                    db.exec(query_direction);
                    if (verbose)
                    {
                        std::cout << "Direction query:" << query_direction << std::endl;
                    }
                }
            }
            transaction.commit();
        }
    }

    void KazSort::Draw(bool verbose)
    {
        std::string cortes_dir = svg_dir_ + "cortes/";
        std::string blocks_dir = svg_dir_ + "blocks/";
        std::string ores_dir = svg_dir_ + "ores/";
        // Draw Cuts | include plan points
        for (const auto &loader_cuts : cuts_by_loader_)
        {
            // Extract loader centroid and outer box
            point_xy min_corner{DBL_MAX, DBL_MAX}, max_corner{-DBL_MAX, -DBL_MAX};
            point_xy loader_centroid;
            // Compute box
            for (const auto &cut : loader_cuts.second)
            {
                float_t x_min{bg::get<0>(min_corner)}, y_min{bg::get<1>(min_corner)};
                float_t x_max{bg::get<0>(max_corner)}, y_max{bg::get<1>(max_corner)};
                auto cut_box = cut.box();
                float_t x_min_cut{bg::get<bg::min_corner, 0>(cut_box)};
                float_t y_min_cut{bg::get<bg::min_corner, 1>(cut_box)};
                float_t x_max_cut{bg::get<bg::max_corner, 0>(cut_box)};
                float_t y_max_cut{bg::get<bg::max_corner, 1>(cut_box)};
                if (x_min_cut < x_min)
                    bg::set<0>(min_corner, x_min_cut);
                if (y_min_cut < y_min)
                    bg::set<1>(min_corner, y_min_cut);
                if (x_max_cut > x_max)
                    bg::set<0>(max_corner, x_max_cut);
                if (y_max_cut > y_max)
                    bg::set<1>(max_corner, y_max_cut);
            }
            bg::set<0>(loader_centroid, (bg::get<0>(max_corner) + bg::get<0>(min_corner)) / 2);
            bg::set<1>(loader_centroid, (bg::get<1>(max_corner) + bg::get<1>(min_corner)) / 2);
            int svg_w = int((bg::get<0>(max_corner) - bg::get<0>(min_corner))) + 100;
            int svg_h = int((bg::get<1>(max_corner) - bg::get<1>(min_corner))) + 100;
            // Draw after apply offset
            std::string svg_file = cortes_dir + loader_cuts.first + ".svg";
            std::string width_height = "width=\"" + std::to_string(svg_w) +
                                       "\" height=\"" + std::to_string(svg_h) + "\"";
            Drawer drawer(svg_file, svg_w, svg_h, width_height);
            if (verbose)
            {
                std::cout << loader_cuts.first << " | centroid ("
                          << loader_centroid.x() << ", "
                          << loader_centroid.y() << ") | box [("
                          << bg::get<0>(min_corner) << ", "
                          << bg::get<1>(min_corner) << "), ("
                          << bg::get<0>(max_corner) << ", "
                          << bg::get<1>(max_corner) << ")] | size ("
                          << svg_w << ", " << svg_h << ")" << std::endl;
            }
            for (int i = 0; i < plan_.size(); i++)
            {
                // apply offset
                point_xy plan_point;
                bg::set<0>(plan_point, (plan_.x()[i] - loader_centroid.x()));
                bg::set<1>(plan_point, (plan_.y()[i] - loader_centroid.y()));
                // draw
                drawer.AddPoint(plan_point);
            }
            for (const auto &cut : loader_cuts.second)
            {
                point_xy cut_moved_centroid;
                MiningCut cut_moved = FixCoordinates(cut, loader_centroid);
                bg::centroid(cut_moved.geometry(), cut_moved_centroid);
                drawer.AddCut(cut_moved, cut.extraction_day());
                drawer.AddNumber(cut.id(), cut_moved_centroid);
                if (verbose)
                {
                    point_xy cut_centroid;
                    bg::centroid(cut.geometry(), cut_centroid);
                    std::cout << "- " << cut
                              << " | centroid ("
                              << cut_centroid.x() << ","
                              << cut_centroid.y() << ") | modified ("
                              << cut_moved_centroid.x() << ","
                              << cut_moved_centroid.y() << ")" << std::endl;
                }
            }
            
        }
        // Draw Blocks
        for (const auto &level_block : blocks_by_level_)
        {
            // Extract loader centroid and outer box
            point_xy min_corner{DBL_MAX, DBL_MAX}, max_corner{-DBL_MAX, -DBL_MAX};
            point_xy level_centroid;
            // Compute box
            for (const auto &block : level_block.second)
            {
                float_t x_min{bg::get<0>(min_corner)}, y_min{bg::get<1>(min_corner)};
                float_t x_max{bg::get<0>(max_corner)}, y_max{bg::get<1>(max_corner)};
                auto block_box = block.geometry();
                float_t x_min_block{bg::get<bg::min_corner, 0>(block_box)};
                float_t y_min_block{bg::get<bg::min_corner, 1>(block_box)};
                float_t x_max_block{bg::get<bg::max_corner, 0>(block_box)};
                float_t y_max_block{bg::get<bg::max_corner, 1>(block_box)};
                if (x_min_block < x_min)
                    bg::set<0>(min_corner, x_min_block);
                if (y_min_block < y_min)
                    bg::set<1>(min_corner, y_min_block);
                if (x_max_block > x_max)
                    bg::set<0>(max_corner, x_max_block);
                if (y_max_block > y_max)
                    bg::set<1>(max_corner, y_max_block);
            }
            bg::set<0>(level_centroid, (bg::get<0>(max_corner) + bg::get<0>(min_corner)) / 2);
            bg::set<1>(level_centroid, (bg::get<1>(max_corner) + bg::get<1>(min_corner)) / 2);
            int svg_w = int((bg::get<0>(max_corner) - bg::get<0>(min_corner))) + 100;
            int svg_h = int((bg::get<1>(max_corner) - bg::get<1>(min_corner))) + 100;
            // Draw after apply offset
            std::string svg_file = blocks_dir + std::to_string(level_block.first) + ".svg";
            std::string width_height = "width=\"" + std::to_string(svg_w) +
                                       "\" height=\"" + std::to_string(svg_h) + "\"";
            Drawer drawer(svg_file, svg_w, svg_h, width_height);
            if (verbose)
            {
                std::cout << level_block.first << " | centroid ("
                          << level_centroid.x() << ", "
                          << level_centroid.y() << ") | box [("
                          << bg::get<0>(min_corner) << ", "
                          << bg::get<1>(min_corner) << "), ("
                          << bg::get<0>(max_corner) << ", "
                          << bg::get<1>(max_corner) << ")] | size ("
                          << svg_w << ", " << svg_h << ")" << std::endl;
            }
            for (const auto &block : level_block.second)
            {
                bool paint = true;
                MiningBlock block_moved = FixCoordinates(block, level_centroid);
                block_moved.set_cu(block.cu());
                drawer.AddBlock(block_moved, paint);
                if (verbose)
                {
                    point_xy block_centroid, block_moved_centroid;
                    bg::centroid(block.geometry(), block_centroid);
                    bg::centroid(block_moved.geometry(), block_moved_centroid);
                    std::cout << "- " << block
                              << " | centroid ("
                              << block_centroid.x() << ","
                              << block_centroid.y() << ") | modified ("
                              << block_moved_centroid.x() << ","
                              << block_moved_centroid.y() << ")" << std::endl;
                }
            }
        }
        // Draw Ores
        for (const auto &level_ores : ores_by_level_)
        {
            // Extract level centroid and outer box
            point_xy min_corner{DBL_MAX, DBL_MAX}, max_corner{-DBL_MAX, -DBL_MAX};
            point_xy level_centroid;
            // Compute box
            for (const auto &ore : level_ores.second)
            {
                float_t x_min{bg::get<0>(min_corner)}, y_min{bg::get<1>(min_corner)};
                float_t x_max{bg::get<0>(max_corner)}, y_max{bg::get<1>(max_corner)};
                auto ore_box = ore.box();
                float_t x_min_ore{bg::get<bg::min_corner, 0>(ore_box)};
                float_t y_min_ore{bg::get<bg::min_corner, 1>(ore_box)};
                float_t x_max_ore{bg::get<bg::max_corner, 0>(ore_box)};
                float_t y_max_ore{bg::get<bg::max_corner, 1>(ore_box)};
                if (x_min_ore < x_min)
                    bg::set<0>(min_corner, x_min_ore);
                if (y_min_ore < y_min)
                    bg::set<1>(min_corner, y_min_ore);
                if (x_max_ore > x_max)
                    bg::set<0>(max_corner, x_max_ore);
                if (y_max_ore > y_max)
                    bg::set<1>(max_corner, y_max_ore);
            }
            bg::set<0>(level_centroid, (bg::get<0>(max_corner) + bg::get<0>(min_corner)) / 2);
            bg::set<1>(level_centroid, (bg::get<1>(max_corner) + bg::get<1>(min_corner)) / 2);
            int svg_w = int((bg::get<0>(max_corner) - bg::get<0>(min_corner))) + 100;
            int svg_h = int((bg::get<1>(max_corner) - bg::get<1>(min_corner))) + 100;
            // Draw after apply offset
            std::string svg_file = ores_dir + std::to_string(level_ores.first) + ".svg";
            std::string width_height = "width=\"" + std::to_string(svg_w) +
                                       "\" height=\"" + std::to_string(svg_h) + "\"";
            Drawer drawer(svg_file, svg_w, svg_h, width_height);
            if (verbose)
            {
                std::cout << level_ores.first << " | centroid ("
                          << level_centroid.x() << ", "
                          << level_centroid.y() << ") | box [("
                          << bg::get<0>(min_corner) << ", "
                          << bg::get<1>(min_corner) << "), ("
                          << bg::get<0>(max_corner) << ", "
                          << bg::get<1>(max_corner) << ")] | size ("
                          << svg_w << ", " << svg_h << ")" << std::endl;
            }
            for (const auto &ore : level_ores.second)
            {
                MiningOre ore_moved = FixCoordinates(ore, level_centroid);
                drawer.AddOre(ore_moved, ore.status());
                if (verbose)
                {
                    point_xy ore_centroid, ore_moved_centroid;
                    bg::centroid(ore.geometry(), ore_centroid);
                    bg::centroid(ore_moved.geometry(), ore_moved_centroid);
                    std::cout << "- " << ore
                              << " | centroid ("
                              << ore_centroid.x() << ","
                              << ore_centroid.y() << ") | modified ("
                              << ore_moved_centroid.x() << ","
                              << ore_moved_centroid.y() << ")" << std::endl;
                }
            }
        }
    }
}