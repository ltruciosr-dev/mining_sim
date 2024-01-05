#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <sstream>
#include "kaz_mining.h"
#include "rapidcsv.h"
#include "SQLiteCpp/SQLiteCpp.h"

namespace astay
{
    KazMining::KazMining(const std::string &db_name) : db_(db_name), transaction_(db_)
    {
    }

    // -- set
    void KazMining::setProjectID(int project_id)
    {
        project_id_ = project_id;
    }

    void KazMining::setProjectDB(const std::string &db_name)
    {
        project_db_ = db_name;
    }

    void KazMining::setModelBlocksDB(const std::string &db_name)
    {
        modelblocks_db_ = db_name;
    }

    void KazMining::setSVGDir(const std::string &svg_dir)
    {
        svg_dir_ = svg_dir;
    }

    // -- read
    void KazMining::ReadCuts(bool verbose)
    {
        SQLite::Database db(project_db_, SQLite::OPEN_READONLY); // project_db
        std::vector<float_t> x_cut, y_cut, z_cut;
        std::vector<int> cut_vert, cut_day, cut_id;
        std::vector<std::string> cut_loader;

        std::string query_cuts = "SELECT X,Y,Z,Vertices,LayerID,Cut_ID,Loader from MiningCuts WHERE Project_ID=" + std::to_string(project_id_) + " ORDER BY IndexPoly ASC, Vertices ASC";
        SQLite::Statement query(db, query_cuts);
        while (query.executeStep())
        {
            float_t x = query.getColumn(0).getDouble();
            float_t y = query.getColumn(1).getDouble();
            float_t z = query.getColumn(2).getDouble();
            int vertice = query.getColumn(3).getInt();
            int day = query.getColumn(4).getInt();
            int id = query.getColumn(5).getInt();
            std::string loader = query.getColumn(6).getString();
            // Populate vectors
            x_cut.push_back(x);
            y_cut.push_back(y);
            z_cut.push_back(z);
            cut_vert.push_back(vertice);
            cut_day.push_back(day);
            cut_id.push_back(id);
            cut_loader.push_back(loader);
        }

        // validate
        if (cut_vert.size() < 3)
        {
            if (verbose)
                std::cout << "[ERROR] ReadCuts | Doesn't have any lines." << std::endl;
            return;
        }
        if (verbose)
            std::cout << "[DEBUG] ReadCuts | Has " << cut_vert.size() << " lines." << std::endl;

        // Offset
        x_offset_ = std::accumulate(x_cut.begin(), x_cut.end(), 0) / x_cut.size();
        y_offset_ = std::accumulate(y_cut.begin(), y_cut.end(), 0) / y_cut.size();
        for (float_t &x : x_cut)
            x -= x_offset_;
        for (float_t &y : y_cut)
            y -= y_offset_;

        // Extract
        int p_vert_n = cut_vert[0] - 1;
        std::string loader;
        MiningCut cut; // a finger has a cut and many ores/blocks
        Finger finger;
        for (size_t i = 0; i < cut_vert.size(); i++)
        {
            int vert_n = cut_vert[i];
            if (vert_n != p_vert_n + 1) // if true, reading new cut
            {
                loader = cut_loader[i - 1];
                if (std::find(loaders_.begin(), loaders_.end(), loader) != loaders_.end())
                {
                    // if loader is in loaders_
                    cut.compute_box();
                    cut.set_level(z_cut[i - 1]);
                    cut.correct();
                    finger.set_cut(cut);
                    finger.set_day(cut_day[i - 1]);
                    finger.set_id(cut_id[i - 1]);
                    fingers_by_loader_[loader].push_back(finger);
                }
                cut.clear();
            }
            point_xy p{x_cut[i], y_cut[i]};
            cut.add(p);
            p_vert_n = vert_n;
        }
        loader = cut_loader.back();
        cut.compute_box();
        cut.set_level(z_cut.back());
        cut.correct();
        finger.set_cut(cut);
        finger.set_day(cut_day.back());
        finger.set_id(cut_id.back());
        fingers_by_loader_[loader].push_back(finger);

        if (verbose)
        {
            std::cout << "Offset: (" << x_offset_ << ", " << y_offset_ << ")" << std::endl;
            std::cout << "Number of Loaders: " << fingers_by_loader_.size() << std::endl;
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -";
            for (const auto &loader_fingers : fingers_by_loader_)
            {
                std::cout << "LOADER [" << loader_fingers.first << "] has (" << loader_fingers.second.size() << ") fingers :" << std::endl;
                int idx = 1;
                for (const auto &finger : loader_fingers.second)
                {
                    std::cout << "Finger ID " << finger.id() << std::endl;
                    std::cout << finger << std::endl;
                    idx++;
                }
                std::cout << std::endl;
            }
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -";
        }
    }

    void KazMining::ReadOres(bool verbose)
    {
        // map (ore_level -> <loader,finger_idx>) | finger_idx is the index of the finger in that lvl.
        std::map<int, std::vector<std::pair<std::string, int>>> fingers_by_level;
        for (auto &loader_fingers : fingers_by_loader_)
        {
            std::string loader = loader_fingers.first;
            int finger_idx = 0;
            for (auto &finger : loader_fingers.second)
            {
                fingers_by_level[finger.cut().level()].emplace_back(loader, finger_idx);
                finger_idx++;
            }
        }

        if (verbose)
        {
            std::cout << "\n-----------------------------------------------" << std::endl;
            std::cout << "ORE HAS " << fingers_by_level.size() << " LEVELS." << std::endl;
            for (const auto &level_fingers : fingers_by_level)
            {
                const int level = level_fingers.first;
                std::cout << "Level " << level << std::endl;
                for (const auto &pair_idx : level_fingers.second)
                {
                    std::string loader = pair_idx.first;
                    int idx = pair_idx.second;
                    std::cout << "- Loader [" << loader << "] | index [" << idx << "] "
                              << " | CutID [" << fingers_by_loader_[loader][idx].id() << "]\n";
                }
            }
        }

        // Read Ores
        SQLite::Database db(modelblocks_db_, SQLite::OPEN_READONLY); // project_db
        for (auto &level_fingers : fingers_by_level)
        {
            // read srg file per level
            int level = level_fingers.first;
            std::vector<float_t> x_ore, y_ore, z_ore;
            std::vector<int> ore_vert;
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
                ore_vert.push_back(vertice);
                ore_name.push_back(name);
            }

            // validate
            if (ore_vert.size() < 3)
            {
                if (verbose)
                    std::cout << "[ERROR] ReadOres | level " << level << " doesn't have any lines." << std::endl;
                continue;
            }
            if (verbose)
            {
                std::cout << "-----------------------------------------------" << std::endl;
                std::cout << "[DEBUG] ReadOres | level " << level << " have " << ore_vert.size() << " lines." << std::endl;
            }

            // substract offset
            for (float_t &x : x_ore)
                x -= x_offset_;
            for (float_t &y : y_ore)
                y -= y_offset_;

            // init extraction
            int p_vert_n = ore_vert[0] - 1;
            MiningOre ore;
            for (size_t i = 0; i < ore_vert.size(); i++)
            {
                int vert_n = ore_vert[i];
                if (vert_n != p_vert_n + 1) // if true, reading new ore
                {
                    ore.compute_box();
                    ore.set_name(ore_name[i - 1]);
                    ore.set_level(z_ore[i - 1]);
                    ore.correct();
                    for (auto &pair_idx : level_fingers.second)
                    {
                        std::string loader = pair_idx.first;
                        int idx = pair_idx.second;
                        if (fingers_by_loader_[loader][idx].ContainsOre(ore))
                            fingers_by_loader_[loader][idx].AddOre(ore);
                    }
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
            for (auto &pair_idx : level_fingers.second)
            {
                std::string loader = pair_idx.first;
                int idx = pair_idx.second;
                if (fingers_by_loader_[loader][idx].ContainsOre(ore))
                    fingers_by_loader_[loader][idx].AddOre(ore);
            }
        }

        if (verbose)
        {
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
            std::cout << "Offset: (" << x_offset_ << ", " << y_offset_ << ")" << std::endl;
            std::cout << "Number of Loaders: " << fingers_by_loader_.size() << std::endl;
            for (const auto &loader_fingers : fingers_by_loader_)
            {
                std::cout << "LOADER [" << loader_fingers.first << "] has (" << loader_fingers.second.size() << ") fingers :" << std::endl;
                int idx = 0;
                for (const auto &finger : loader_fingers.second)
                {
                    std::cout << "Finger " << idx << " | Cut ID " << finger.id() << " | Level " << finger.cut().level() << std::endl;
                    std::cout << finger << std::endl;
                    idx++;
                }
                std::cout << std::endl;
            }
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -";
        }
    }

    void KazMining::ReadBlocks(bool verbose)
    {
        // Extract model
        ModelBlock model_block(modelblocks_db_);
        model_block.set_offset(x_offset_, y_offset_);
        model_block.read_setup();
        if (verbose)
        {
            std::cout << "--------------------------------------" << std::endl;
            std::cout << "Model Block: " << std::endl;
            std::cout << model_block << std::endl;
        }
        for (auto &loader_fingers : fingers_by_loader_)
        {
            for (auto &finger : loader_fingers.second)
            {
                std::vector<MiningBlock> blocks = model_block.get_blocks(finger);
                finger.AddBlocks(blocks);
            }
        }

        if (verbose)
        {
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - " << std::endl;
            std::cout << "Offset: (" << x_offset_ << ", " << y_offset_ << ")" << std::endl;
            std::cout << "Number of Loaders: " << fingers_by_loader_.size() << std::endl;
            for (const auto &loader_fingers : fingers_by_loader_)
            {
                std::string loader = loader_fingers.first;
                std::cout << "LOADER [" << loader << "] has (" << loader_fingers.second.size() << ") fingers :" << std::endl;
                int idx = 0;
                for (const auto &finger : loader_fingers.second)
                {
                    std::cout << "Finger " << idx << " | Cut ID " << finger.id() << " | Level " << finger.cut().level() << std::endl;
                    std::cout << finger << std::endl;
                    idx++;
                }
                std::cout << std::endl;
            }
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -";
        }
    }

    void KazMining::ComputeSlices(float_t slice_dist, bool verbose)
    {
        SQLite::Database db(project_db_, SQLite::OPEN_READONLY); // project_db
        for (auto &loader_fingers : fingers_by_loader_)
        {
            // Avoid compute if loader not in loaders_
            std::string loader = loader_fingers.first;
            if (std::find(loaders_.begin(), loaders_.end(), loader) != loaders_.end())
            {
                if (verbose)
                    std::cout << "Valid Loader > [" << loader << "]" << std::endl;
            }
            else
            {
                if (verbose)
                    std::cout << "Invalid Loader > [" << loader << "]" << std::endl;
                continue;
            }

            // Get direction vector and store by finger_id of certain loader
            std::unordered_map<int, point_xy> direction_by_id;
            std::string query_direction = "SELECT Cut_ID,Dir_X,Dir_Y from MiningDirection WHERE Project_ID=" + std::to_string(project_id_) +
                                          " AND Shovel='" + loader + "'";
            SQLite::Statement query(db, query_direction);
            while (query.executeStep())
            {
                // -- get
                int id = query.getColumn(0).getInt();
                float_t dir_x = query.getColumn(1).getDouble();
                float_t dir_y = query.getColumn(2).getDouble();
                // -- process
                point_xy dir_vec{dir_x, dir_y};
                float l2 = sqrt(bg::get<0>(dir_vec) * bg::get<0>(dir_vec) +
                                bg::get<1>(dir_vec) * bg::get<1>(dir_vec));
                if (l2 > 0.1)
                {
                    bg::divide_value(dir_vec, l2);
                    direction_by_id[id] = dir_vec;
                }
            }

            // Compute slices & limits
            int N = loader_fingers.second.size();
            point_xy ref_point{0, 0};
            std::vector<int> values;
            std::vector<float_t> v_limits;
            for (auto &finger : loader_fingers.second)
            {
                int id = finger.id();
                point_xy dir_vec;
                float_t limit;
                if (direction_by_id.count(id) > 0) // direction is set
                {
                    dir_vec = direction_by_id[id];
                    // -- slice
                    finger.ComputeSlices(slice_dist, dir_vec);
                    auto limits = finger.ComputeLimitsByDirection(dir_vec, ref_point);
                    limit = limits[0];
                }
                else // direction is not set
                {
                    bg::set<0>(dir_vec, 0.0);
                    bg::set<1>(dir_vec, 0.0);
                    limit = 10E5;
                }
                // -- limits
                v_limits.push_back(limit);
                values.push_back(finger.day() * N);
                // -- verbose
                if (verbose)
                {
                    std::cout << "> [" << loader
                              << "] | Finger: ID(" << id
                              << "), Day (" << finger.day()
                              << "), Limit (" << limit
                              << ") | Direction [" << bg::get<0>(dir_vec)
                              << ", " << bg::get<1>(dir_vec) << "]" << std::endl;
                }
            }

            std::vector<int> v_order(N);
            std::iota(v_order.begin(), v_order.end(), 0);
            std::sort(v_order.begin(), v_order.end(), [&](int i, int j)
                      { return v_limits[i] < v_limits[j]; });

            // Sort mining order
            int idx = 0;
            for (int order : v_order)
                values[order] += idx++;
            std::iota(v_order.begin(), v_order.end(), 0);
            std::sort(v_order.begin(), v_order.end(), [&](int i, int j)
                      { return values[i] < values[j]; });
            fingers_order_by_loader_[loader] = std::move(v_order);
        }
        if (verbose)
        {
            std::cout << "Number of Loaders: " << fingers_by_loader_.size() << std::endl;
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -";
            for (const auto &loader_fingers : fingers_by_loader_)
            {
                std::string loader = loader_fingers.first;
                std::cout << "LOADER [" << loader << "] has ("
                          << loader_fingers.second.size() << ") fingers :" << std::endl;
                int idx = 0;
                for (int finger_order : fingers_order_by_loader_[loader])
                {
                    std::cout << "Order [" << idx << "] | Cut ID "
                              << loader_fingers.second[finger_order].id() << std::endl;
                    idx++;
                }
            }
            std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -";
        }
    }

    void KazMining::ComputeBlockIntersections(const polygon_t &ore_inter, int index, SQLite::Database &db)
    {
        t_area_ = 0;
        t_cu_ = 0;
        t_ton_ = 0;
        const auto ore_inter_box = Polygon::GetBoxFromPolygon(ore_inter);
        for (const auto &block : fingers_by_loader_[loader_][index].blocks())
        {
            if (bg::intersects(ore_inter_box, block.geometry()))
            {
                std::deque<polygon_t> block_inters;
                bg::intersection(block.geometry(), ore_inter, block_inters);
                for (const auto &block_inter : block_inters)
                {
                    float_t area = bg::area(block_inter);
                    float_t cu = block.cu() * area;
                    float_t ton = block.density() * area * 10;
                    t_area_ += area;
                    t_cu_ += cu;
                    t_ton_ += ton;
                }
            }
        }
        t_cu_ /= t_area_;
        std::stringstream buffer;
        buffer << "INSERT INTO Autoslicer (Project_ID, Loader,Cut_ID,Slice,OreType,Ore_X,Ore_Y,Ore_Z,Area,Cu,Ton) VALUES ("
               << project_id_ << ", '"
               << loader_ << "', "
               << finger_id_ << ", "
               << num_slice_ << ", '"
               << ore_name_ << "', "
               << ore_position_[0] << ", " << ore_position_[1] << ", " << ore_position_[2] << ", "
               << t_area_ << ", " << t_cu_ << ", " << t_ton_ << ")";

        db.exec(buffer.str());
    }

    void KazMining::ComputeOreIntersections(const polygon_t &cut_inter, int index, SQLite::Database &db)
    {
        const auto cut_inter_box = Polygon::GetBoxFromPolygon(cut_inter);
        for (const auto &ore : fingers_by_loader_[loader_][index].ores())
        {
            if (bg::intersects(cut_inter_box, ore.box()))
            {
                std::deque<polygon_t> ore_inters;
                bg::intersection(ore.geometry(), cut_inter, ore_inters);
                if (!ore_inters.empty())
                {
                    ore_name_ = ore.name();
                    for (const auto &ore_inter : ore_inters)
                    {
                        point_xy ore_centroid;
                        boost::geometry::centroid(ore_inter, ore_centroid);
                        ore_position_[0] = bg::get<0>(ore_centroid) + x_offset_;
                        ore_position_[1] = bg::get<1>(ore_centroid) + y_offset_;
                        ore_position_[2] = ore.level();
                        ComputeBlockIntersections(ore_inter, index, db);
                    }
                }
            }
        }
    }

    void KazMining::ComputeFingerIntersections(bool verbose)
    {
        SQLite::Database db(project_db_, SQLite::OPEN_READWRITE);
        SQLite::Transaction transaction(db);
        // -- remove autoslicer
        std::string query_remove = "DELETE FROM Autoslicer WHERE Project_ID=" + std::to_string(project_id_);
        db.exec(query_remove);
        for (auto &loader_order : fingers_order_by_loader_)
        {
            loader_ = loader_order.first; // set loader in iteration
            for (const int index : loader_order.second)
            {
                const auto &cut = fingers_by_loader_[loader_][index].cut();
                finger_id_ = fingers_by_loader_[loader_][index].id(); // set cut in iteration
                int slice_idx = 0;
                for (const auto slice : fingers_by_loader_[loader_][index].slices())
                {
                    num_slice_ = slice_idx++;
                    std::deque<polygon_t> cut_inters;
                    bg::intersection(cut.geometry(), slice.geometry(), cut_inters);
                    if (!cut_inters.empty())
                    {
                        for (const auto &cut_inter : cut_inters)
                        {
                            ComputeOreIntersections(cut_inter, index, db);
                        }
                    }
                }
            }
        }
        transaction.commit();
    }

    void KazMining::Draw(bool verbose)
    {
        // Draw Fingers | autoslicer/
        for (auto &loader_fingers : fingers_by_loader_)
        {
            // Avoid compute if loader not in loaders_
            std::string loader = loader_fingers.first;
            if (std::find(loaders_.begin(), loaders_.end(), loader) != loaders_.end())
            {
                if (verbose)
                    std::cout << "Valid Loader > [" << loader << "]" << std::endl;
            }
            else
            {
                if (verbose)
                    std::cout << "Invalid Loader > [" << loader << "]" << std::endl;
                continue;
            }

            // Extract loader centroid and outer box
            point_xy min_corner{DBL_MAX, DBL_MAX}, max_corner{-DBL_MAX, -DBL_MAX};
            point_xy loader_centroid;

            // Compute box
            for (auto &finger : loader_fingers.second)
            {
                const auto &cut = finger.cut();
                float_t x_min{bg::get<0>(min_corner)}, y_min{bg::get<1>(min_corner)};
                float_t x_max{bg::get<0>(max_corner)}, y_max{bg::get<1>(max_corner)};
                const auto &cut_box = cut.box();
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
            int svg_w = int((bg::get<0>(max_corner) - bg::get<0>(min_corner))) + 150;
            int svg_h = int((bg::get<1>(max_corner) - bg::get<1>(min_corner))) + 150;

            // Draw after apply offset
            std::string svg_file = svg_dir_ + "autoslicer/" + loader + ".svg";
            std::string width_height = "width=\"" + std::to_string(svg_w) +
                                       "\" height=\"" + std::to_string(svg_h) + "\"";
            Drawer drawer(svg_file, svg_w, svg_h, width_height);
            if (verbose)
            {
                std::cout << loader << " | centroid ("
                          << loader_centroid.x() << ", "
                          << loader_centroid.y() << ") | box [("
                          << bg::get<0>(min_corner) << ", "
                          << bg::get<1>(min_corner) << "), ("
                          << bg::get<0>(max_corner) << ", "
                          << bg::get<1>(max_corner) << ")] | size ("
                          << svg_w << ", " << svg_h << ")" << std::endl;
            }
            for (auto &finger : loader_fingers.second)
            {
                const auto &cut = finger.cut();
                point_xy cut_moved_centroid;
                MiningCut cut_moved = FixCoordinates(cut, loader_centroid);
                bg::centroid(cut_moved.geometry(), cut_moved_centroid);
                // drawer.AddCut(cut_moved, cut.extraction_day());
                // drawer.AddNumber(cut.id(), cut_moved_centroid);
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
    }
}